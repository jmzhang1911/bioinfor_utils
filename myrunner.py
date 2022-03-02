#!/usr/bin/env python
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
import functools
import datetime
import logging
import random
import time

FORMAT = '%(asctime)s %(threadName)s=> %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class MyPath:
    @staticmethod
    def mkdir(*path: str):
        """mkdir recursive-dir"""
        for p in list(path):
            if not Path(p).exists():
                Path(p).mkdir(parents=True, exist_ok=True)
            logging.info('{} has been created'.format(p))


class MyRunner:
    """
    1) Using ThreadPoolExecutor to run shell commands
    There are two ways to do it:
        - @MyRunner.cmd_wrapper(): the wrapper-ed function must return a list contains shell commands
        - MyRunner.runner(cmd_list): just use it as a function

    Notes:
        - using absolute path when using qsub
        - mkdir a dir named __work_sh for saving cmds
        - mkdir a dir named ____qsub_sh for saving qsub cmds
        - threads: running all cmds at same time(default, using threads_num to set threads num)
        - test: for testing the programs not actually running a linux cmd just print the cmds
        - cmd_name: set cmds name default(func.__name__)

    2) count the running time of function
    ...@count_running_time
    ...def foo():
    ...    return

    """
    run_qsub = '/share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh'

    @staticmethod
    def check_file(*file_list):
        logging.info('checking the input files existed or not ...')
        for i in file_list:
            if not Path(i).exists():
                raise '{} is not existed'.format(i)
        logging.info('input files checked -> done')

    @staticmethod
    def run_shell(cmd, test=False):
        logging.info('subprocess is running cmd={}'.format(cmd))
        if not test:
            p = subprocess.Popen(cmd, shell=True)
            ret = p.wait()
            if ret != 0:
                raise Exception('cmd=`{}` failed:\n{}'.format(cmd, ret))
        else:
            time.sleep(random.randint(4, 10))

    @staticmethod
    def _run_cmd(func_name, cmd_list):
        if not isinstance(cmd_list, list):
            return

        MyPath.mkdir('cmds_sh')
        cmd_file = 'cmds_sh/{}.sh'.format(func_name)
        with open(cmd_file, 'w') as f:
            for cmd in cmd_list:
                f.write(cmd + '\n')

    # TODO:此方法基于perl+shell实现，未来尝试直接使用python实现
    def _run_qsub(self, func_name, cmd_list, queue):
        MyPath.mkdir('qsub_sh')
        func_dir = './qsub_sh/{}'.format(func_name)
        MyPath.mkdir(func_dir)
        tmp_list = []
        for num, cmd in enumerate(cmd_list):
            cmd_file = Path(func_dir) / (func_name + '_' + str(num) + '.sh')
            with open(cmd_file, 'w') as f:
                f.write(cmd + '\n')

            cmd = 'bash {} --reqsub --independent {} --queue {}'.format(self.run_qsub, cmd_file, queue)
            tmp_list.append(cmd)

        return cmd_list

    @classmethod
    def cmd_wrapper(cls, cmd_name='', test=False, qsub=False, queue='low.q', threads_num=''):
        def _wrapper(fn):
            @functools.wraps(fn)
            def __wrapper(*args, **kwargs):
                cmd_list = fn(*args, **kwargs)
                save_name = fn.__name__ if not cmd_name else cmd_name
                cls.runner(cmd_list,
                           func_name=save_name,
                           test=test,
                           qsub=qsub,
                           queue=queue,
                           threads_num=threads_num)

            return __wrapper

        return _wrapper

    @classmethod
    def runner(cls, cmd_list, func_name='', qsub=False, threads_num='', test=False, queue='low.q'):
        cmd_list = cmd_list[0] if isinstance(cmd_list, tuple) else cmd_list
        try:
            threads_num = len(cmd_list) if not threads_num else threads_num
        except TypeError:
            threads_num = 1

        if not cmd_list:
            logging.info('myrunner|> Exception: there are no cmds')

        if qsub:
            cmd_list = cls._run_qsub(func_name, cmd_list, queue)
        else:
            cls._run_cmd(func_name, cmd_list)

        with ThreadPoolExecutor(max_workers=threads_num) as executor_cmd:
            run_shell = functools.partial(cls.run_shell, test=test)
            task_list = []
            for cmd in cmd_list:
                task = executor_cmd.submit(run_shell, cmd)
                task_list.append(task)

            for future in as_completed(task_list):
                future.result()

    @classmethod
    def count_running_time(cls, fn):
        @functools.wraps(fn)
        def _wrapper(*args, **kwargs):
            start = datetime.datetime.now()
            res = fn(*args, **kwargs)
            delta = (datetime.datetime.now() - start).total_seconds()
            logging.info('%%%%%%%%-- {}.{} took {:.2} min --%%%%%%%%'.format(fn.__module__, fn.__name__, delta / 60))
            return res

        return _wrapper
