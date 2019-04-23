import time
import random
import hashlib
from numpy import base_repr
from threading import Thread, Lock
from multiprocessing import Process, Queue, current_process, freeze_support

class Job:
    def __init__(self, func, args):
        self.id = 0
        self.func = func
        self.args = args
        self.queue = Queue()
        self.process = Process(target=self.calculate, args=())
        self.started = False
        self.finished = False
        self.queuetime = time.time()
        self.starttime = 0
        self.finishtime = 0
        self.error = False

    def calculate(self):
        try:
            res = self.func(*self.args)
            self.queue.put(res)
            self.queue.task_done()
        except:
            self.error = True

class JobManager:
    def __init__(self, maxprocesses):
        self.mutex = Lock()
        self.maxprocesses = maxprocesses
        self.activeprocesses = 0
        self.jobcount = 0
        self.jobs = {}
        self.cancelqueue = Queue()

    def addjob(self, func, args):
        self.mutex.acquire()
        jobid = 0
        try:
            self.jobcount += 1
            job = Job(func,args)
            s = "%d" % self.jobcount
            job.id = hashlib.sha256(s.encode('utf-8')).hexdigest()
            self.jobs[job.id] = job
            print("JOB %s" % job.id)
            jobid = job.id
        finally:
            self.mutex.release()
        self.update()
        return jobid

    def canceljob(self, jobid):
        self.mutex.acquire()
        try:
            self.cancelqueue.put(jobid)
        finally:
            self.mutex.release()


    def update(self):
        self.mutex.acquire()
        try:
            print("Jobs: %d" % self.jobcount)
            while not self.cancelqueue.empty():
                jobid = self.cancelqueue.get()
                if jobid in self.jobs:
                    job = self.jobs[jobid]
                    if job.started:
                        job.process.terminate()
                    del self.jobs[jobid]

            for jobid in self.jobs.keys():
                job = self.jobs[jobid]
                if not job.started and self.activeprocesses < self.maxprocesses:
                    job.process.start()
                    job.started = True
                    job.starttime = time.time()
                    self.activeprocesses += 1
                if job.started and not job.finished and not job.process.is_alive():
                    job.finished = True
                    job.finishtime = time.time()
                    self.activeprocesses -= 1
                    print("GETING %s %s" % (job.id, job.process.exitcode))
                    try:
                        print("RES %s" % job.queue.get_nowait())
                    except:
                        print("EXCEPTION IN %s" % job.id)
                        job.process.terminate()
                print("ID %s Started %s, Finished %s" % (job.id, job.started, job.finished))
        finally:
            self.mutex.release()

    def loop(self):
        while True:
            self.update()
            time.sleep(1.0)
