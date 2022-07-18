from bioinfor_tools.cmd_runner import CmdRunner

@CmdRunner.cmd_wrapper()
def foo():
    return ['echo {}'.format(_) for _ in range(10)]

foo()

CmdRunner.cmd(['echo 2'])