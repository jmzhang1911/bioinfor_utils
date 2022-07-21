import logging
from pathlib import Path

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')


class BioBasic:
    """
    - basic functions and property
        - read cfg
        - mkdir dir
        - module property
    """
    # 主路径位置
    BASE_DIR = Path(__file__).absolute().parent.parent
    # 日志路径
    LOGGING_SUMMARY = BASE_DIR / 'logging_summary.txt'

    def __init__(self, module):
        # 模块所在的绝对父目录
        self._module_dir = Path(module).absolute().parent
        # 模块的名称
        self._module_name = Path(module).name

    @staticmethod
    def mkdir(*path: str):
        """mkdir recursive-dir"""
        for p in list(path):
            if not Path(p).exists():
                Path(p).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def read_cfg(config):
        with open(config, 'r', encoding='utf-8') as f:
            config_dict = {}
            for line in f:
                if line.strip().startswith('#') or line.strip() == '':
                    continue
                k, v = line.strip().split()
                config_dict[k] = str(v).strip()
        return config_dict
