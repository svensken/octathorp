import pp
from multiprocessing import Process

def miin(a):
    return min(a)


if __name__ == '__main__':
    job_server = pp.Server()
    print job_server.get_ncpus()
    
    a = [1,2,3,4,5,6,7,8,9,10]

    job1 = job_server.submit(miin, (a,),(),())
    job2 = job_server.submit(miin, (a,),(),())
    job3 = job_server.submit(miin, (a,),(),())

    j1 = job1()
    j2 = job2()
    j3 = job3()
