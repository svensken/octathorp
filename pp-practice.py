import pp
import time

def miin(a):
    time.sleep(20)
    return min(a)

def maax(a):
    time.sleep(3)
    return max(a)


if __name__ == '__main__':
    job_server = pp.Server()
    print job_server.get_ncpus()
    
    a = [1,2,3,4,5,6,7,8,9,10]

    job1 = job_server.submit(miin, (a,), (), ('time',) )
    job2 = job_server.submit(maax, (a,), (), ('time',) )

    j1 = job1()
    j2 = job2()
    print j2
    print j1
