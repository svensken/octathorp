from progressbar import ProgressBar
import time

prog1 = ProgressBar().start()
prog2 = ProgressBar().start()

for i in range(100):
    for h in range(100):
        time.sleep(.01)
        print 'fdsa'
        prog2.update(h+1)
    prog1.update(i+1)
