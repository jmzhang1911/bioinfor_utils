from myrunner import *
import threading
import time

# @MyRunner.cmd_wrapper(test=True, threads_num=5)
# def foo():
#     time.sleep(2)
#     return ['sleep {}'.format(str(i)) for i in range(10)]
#
#
# @MyRunner.cmd_wrapper(test=True)
# def foo2():
#     return ['kill {}'.format(str(i)) for i in range(10)]
#
#
# foo()
# foo2()


print('a{0} b{0} c{0} d{0}'.format(('abc ',)))
# print(*('abc',) * 4)