from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
import functools
import datetime
import getpass
import logging
import time

FORMAT = '%(asctime)s %(threadName)s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class MyPath:
    @staticmethod
    def mkdir(*path: str):
        """mkdir recursive-dir"""
        for p in list(path):
            if not Path(p).exists():
                Path(p).mkdir(parents=True, exist_ok=True)


class CmdRunner:
    """
    1) basic
    Note: cmd_container: ([cmd1, cmd2, *], *) or [cmd1, cmd2, *]
    There are two ways to do it:
        - @CmdRunner.cmd_wrapper(): the wrapper-ed function must return a cmd_container
        - CmdRunner.runner(cmd_container): just use it as a function

        example 01:

        @CmdRunner.cmd_wrapper()
        def hello_world():
            cmd_list = ['sleep 5 && echo done' for _ in range(5)]
            return cmd_list

        example 02:
        CmdRunner.runner(['echo hello world'])

    2) Using PBS: add parameter use_qsub=True
    3) outputs:
        - cmds_sh:
            - func_name-time.sh # contains cmds
        - cmds_qsub:
            - func_name-time
                - cmd_1.pbs # pbs script
                - cmd_1.xxx # means failure
                - cmd_1.output # standard output
                - cmd_1.error # standard error output
                - cmd_1.done # means success
    """

    user_name = getpass.getuser()
    qsub_nodes = 'hpc12-l-0-0'

    # qsub的一些属性
    _qsub_host = 'ssh {}@{} "cd {} && '.format(user_name, qsub_nodes, Path().absolute())
    _qsub_re_submit_times = 3  # 重投次数
    _qsub_check_time = 60  # 多久检查一次任务执行状况，单位：秒

    @staticmethod
    def _resolve_cmd(cmd_container):
        """
        解析cmd -> ([cmd1, cmd2, *], *) or [cmd1, cmd2, *]
        """
        if not cmd_container:
            raise 'None cmd, shut down!!!'

        if isinstance(cmd_container, list):
            return cmd_container, None

        elif isinstance(cmd_container, tuple) and isinstance(cmd_container[0], list):
            return cmd_container[0], cmd_container[1:]

        else:
            raise Exception('please input ([cmd1, cmd2, *], *) or [cmd1, cmd2, *]')

    @classmethod
    def _run_shell(cls, cmd, use_qsub: bool, file_exist_for_skipping=''):
        """
        - 执行命令，可选普通shell模式或qsub模式
        - qsub模式，通过qstat判断任务是否正在进行，运行完成后判断是否touch.done，若失败会重新运行三次
        """
        if not use_qsub:
            logging.info('-> running {}'.format(cmd))
            p = subprocess.Popen(cmd, shell=True)
            ret = p.wait()

            if ret != 0:
                raise Exception('-> `{}` failed:\n{}'.format(cmd, ret))

        else:
            cmd_qsub, touch_done, cmd_shell = cmd
            num = 1

            while num <= cls._qsub_re_submit_times:
                logging.info('-> qsub is running {}'.format(cmd_shell))
                p = subprocess.Popen(cmd_qsub, shell=True, stdout=subprocess.PIPE)
                p.wait()
                job_number = p.stdout.read().decode().split()[2]

                while True:
                    time.sleep(cls._qsub_check_time)
                    p1 = subprocess.Popen('qstat', shell=True, stdout=subprocess.PIPE)
                    msg = p1.stdout.read().decode()

                    if job_number not in [_.split()[0] for _ in msg.split('\n')[2:-1]]:

                        if Path(touch_done.replace('&& touch ', '')).exists():
                            num = 666
                            logging.info('cmd={} done! Good for U'.format(cmd_shell))
                            break
                        else:
                            num += 1
                            logging.info('cmd={} shut down!!!trying {} times ...'.format(cmd_shell, num))
                            if num > cls._qsub_re_submit_times:
                                with open(cmd_shell + '.xxx', 'w') as f:
                                    f.write('sorry for you, {}!!!\n'.format(cls.user_name))
                                raise Exception('cmd={} failed!!!tried {} times ...'.format(cmd_shell, num))

    @classmethod
    def cmd(cls, cmd_container, n_jobs: int = '', fn_name='', module_name='', use_qsub=False):

        logging.info(
            '>>>>-- Dear {}, {}.{} is ready to go. Godspeed!!! --<<<<'.format(cls.user_name, module_name, fn_name))
        # 解析出cmd_list，函数的返回值等
        cmd_list, return_value = cls._resolve_cmd(cmd_container)

        n_jobs = int(n_jobs) if n_jobs else len(cmd_list)
        now = str(datetime.datetime.now().replace(microsecond=0)).replace(' ', '_').replace(':', '-')

        if not use_qsub:
            if fn_name:
                MyPath.mkdir('cmds_sh')
                cmd_file = 'cmds_sh/{}_{}.sh'.format(fn_name, now)

                with open(cmd_file, 'w') as f:
                    for cmd in cmd_list:
                        f.write(cmd + '\n')
        else:
            fn_name = fn_name if fn_name else 'tmp_qsub_cmd'
            fn_path = 'cmds_qsub/{}_{}'.format(fn_name, now)
            MyPath.mkdir(fn_path)
            cmd_qsub_list = []

            for num, cmd in enumerate(cmd_list, start=1):
                cmd_shell = '{}/{}_{}.pbs'.format(fn_path, fn_name, num)
                cmd_shell_o = cmd_shell + '.output'
                cmd_shell_e = cmd_shell + '.error'

                with open(cmd_shell, 'w') as f:
                    touch_done = '&& touch {}.done'.format(cmd_shell)
                    cmd += touch_done
                    f.write(cmd + '\n')

                # 解决在计算节点投递任务的情况

                cmd_qsub = '{} qsub -V -cwd -o {} -e {} {}"'.format(cls._qsub_host, cmd_shell_o, cmd_shell_e,
                                                                    cmd_shell)
                cmd_qsub_list.append((cmd_qsub, touch_done, cmd_shell))

            cmd_list = cmd_qsub_list

        with ThreadPoolExecutor(max_workers=n_jobs) as executor_cmd:
            run_shell = functools.partial(cls._run_shell, use_qsub=use_qsub)
            task_list = []
            for cmd in cmd_list:
                task = executor_cmd.submit(run_shell, cmd)
                task_list.append(task)

            for future in as_completed(task_list):
                future.result()

        return return_value

    # def __call__(self, n_jobs: int = '', *args, use_qsub=False, **kwargs):
    #     start = datetime.datetime.now()
    #     cmd_container = self._resolve_cmd(self.fn(*args, **kwargs))
    #     return_value = self.cmd(cmd_container=cmd_container,
    #                             n_jobs=n_jobs,
    #                             use_qsub=use_qsub,
    #                             fn_name=self.fn.__name__,
    #                             module_name=self.fn.__module__)
    #
    #     delta = (datetime.datetime.now() - start).total_seconds()
    #     logging.info('%%%%-- {}.{} took {:.1} min --%%%%'.format(self.fn.__module__, self.fn.__name__, delta / 60))
    #     return return_value

    @classmethod
    def cmd_wrapper(cls, n_jobs: int = '', use_qsub=False):

        def wrapper(fn):
            @functools.wraps(fn)
            def __wrapper(*args, **kwargs):
                start = datetime.datetime.now()
                cmd_container = cls._resolve_cmd(fn(*args, **kwargs))
                return_value = cls.cmd(cmd_container=cmd_container,
                                       n_jobs=n_jobs,
                                       use_qsub=use_qsub,
                                       fn_name=fn.__name__,
                                       module_name=fn.__module__)

                delta = (datetime.datetime.now() - start).total_seconds()
                logging.info('%%%%-- {}.{} took {:.1} min --%%%%'.format(fn.__module__, fn.__name__, delta / 60))

                return return_value

            return __wrapper

        return wrapper


if __name__ == '__main__':
    @CmdRunner.cmd_wrapper()
    def foo(a=10):
        cmd_list = ['sleep 3 && echo done!' for _ in range(20)]
        return cmd_list, a


    import random


    @CmdRunner.cmd_wrapper(use_qsub=True, n_jobs=5)
    def foo2(a='hello'):
        cmd_list = ['sleep {} && echo done!'.format(random.randint(15, 30)) for _ in range(10)]
        cmd_list.append('xxxxx ')
        return cmd_list, a


    @CmdRunner.cmd_wrapper(use_qsub=True)
    def hello_world(a='hello'):
        cmd_list = ['sleep {} && echo done!'.format(random.randint(15, 30)) for _ in range(10)]
        cmd_list.append('xxxxx ')
        return cmd_list, a

    # @CmdRunner(n_jobs=1)  # mdRunner()
    # def hello_world2(a='hello'):
    #     cmd_list = ['sleep {} && echo done!'.format(random.randint(15, 30)) for _ in range(10)]
    #     cmd_list.append('xxxxx ')
    #     return cmd_list, a
