
import threading
import time

ls = []
BER = [0]
BLER = [0]
def count(start, stop):
    print("hey")
    for i in range(start, stop+1):
        ls.append(i)
        BER[0] = BER[0] + i
        BLER[0] = BLER[0] + 1

def timer():
    print("heyo")
    time.sleep(10)
    print("ok bye")
x = threading.Thread(target=count, args=(0, 5))
y = threading.Thread(target=timer, args=())

x.start()
y.start()

x.join()
y.join()
print(ls)
print(BER)
print(BLER)
print(sum(ls))
"""
def foo(bar, baz):
  print('hello {0}'.format(bar))
  return 'foo' + baz

from multiprocessing.pool import ThreadPool
pool = ThreadPool(processes=2)

async_result = pool.apply_async(foo, ('world', 'foo')) # tuple of args for foo
asi = pool.apply_async(foo, ('meta', 'verse'))

# do some other stuff in the main process

return_val = async_result.get()  # get the return value from your function.
ret = asi.get()
"""
