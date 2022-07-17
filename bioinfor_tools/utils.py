from pathlib import Path
import datetime


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
