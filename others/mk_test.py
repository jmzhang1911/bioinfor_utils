from pathlib import Path
import logging
from bioinfor_tools.utils import BioUtils


class Tmp(BioUtils):

    def __init__(self):
        super().__init__(module=__file__)
        logging.info('xi')

    pass



t = Tmp()
t.copy_readme()
t.make_summary()
print(t.__dict__)
print(t.__class__.__dict__)

for i in ['asdf asdf']:
    logging.info(i)
