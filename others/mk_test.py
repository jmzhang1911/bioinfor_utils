from bioinfor_tools.cmd_runner import CmdRunner


@CmdRunner.cmd_wrapper(n_jobs=5)
def foo():
    return ['sleep {}'.format(_) for _ in range(2)]

foo()

print(foo.__module__.__class__.__file__)