from bioinfor_tools.cmd_runner import CmdRunner
from bioinfor_tools.bio_basic import BioBasic
from pathlib import Path

import shutil
import logging
import datetime


class BioUtils(BioBasic):
    """
    - tools for pipeline
        - copy readme
        - write log
        - mv results
        - etc..
    """
    # {module_name : readme}
    README = {'run_GSEA_GSVA.py': ['BMK_GSEA结果说明文档202206.pdf', 'BMK_GSVA结果说明文档202206.pdf'],
              'run_infercnv.py': ['BMK_inferCNV结果说明文档202206.pdf'],
              'run_scenic.py': ['BMK_pySCIENE结果说明文档202206.pdf'],
              'run_cellphonedb.py': ['BMK_细胞通讯分析结果说明文档202206.pdf']}

    def __init__(self, module):
        super().__init__(module)
        self._status = 0

    def copy_readme(self, output='.'):
        readme_dict = {k: '{}/README_FILE/{}'.format(self.BASE_DIR, v) for k, v in self.README.items()}
        src = readme_dict.get(self._module_name)
        logging.info(src)
        if src and Path(src).exists() and Path(output).exists():
            for file in src:
                shutil.copy(file, str(output))
        else:
            logging.info('Sorry, readme not found or output not found ...')

    def make_summary(self):
        status = 'doing' if self._status == 0 else 'done'

        with open(self.LOGGING_SUMMARY, 'a') as f:
            time_now = datetime.datetime.now().replace(microsecond=0)
            info_log = 'datetime:{}\tpwd:{}\tscript_id:{}\tstatus:{}\n'. \
                format(time_now, Path().absolute(), self._module_name, status)
            f.write(info_log)

            self._status += 1

    @CmdRunner.cmd_wrapper(n_jobs=5)
    def use_config_to_copy(self, copy_config):
        return list()


class MyPath:
    @staticmethod
    def mkdir(*path: str):
        """mkdir recursive-dir"""
        for p in list(path):
            if not Path(p).exists():
                Path(p).mkdir(parents=True, exist_ok=True)


def make_summary(script_path, status='doing'):
    logging_summary = Path(script_path).parent.parent / 'logging_summary.txt'
    with open(logging_summary, 'a') as f:
        info_log = 'datetime:{}\tpwd:{}\tscript_id:{}\tstatus:{}\n'. \
            format(datetime.datetime.now().replace(microsecond=0), Path().absolute(), Path(script_path).name, status)
        f.write(info_log)


def read_cfg(config):
    with open(config, 'r', encoding='utf-8') as f:
        config_dict = {}
        for line in f:
            if line.strip().startswith('#') or line.strip() == '':
                continue
            k, v = line.strip().split()
            config_dict[k] = str(v).strip()
    return config_dict
