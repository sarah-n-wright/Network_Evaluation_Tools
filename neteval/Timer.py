from datetime import datetime

class Timer:
    def __init__(self):
        self.start_times = {}
        self.finish_times = {}
        self.elapsed_times = {}
        self.tasks = []
        self.current_task_stack = []
        self.indents = {}
        
    def start(self, taskstr):
        if taskstr in self.start_times.keys():
            taskstr=taskstr + "1"
            i=1
            while taskstr in self.start_times.keys():
                i += 1
                taskstr = taskstr[0:-1] + str(i) 
        self.current_task_stack.append(taskstr)
        self.indents[taskstr] = len(self.current_task_stack) - 1
        self.tasks.append(taskstr)
        self.start_times[taskstr] = datetime.now()
        
    def end(self, taskstr):
        if taskstr in self.finish_times:
            matching_tasks = [task for task in self.start_times.keys() if taskstr in task]
            taskstr = matching_tasks[-1]
        self.current_task_stack.remove(taskstr)
        self.finish_times[taskstr] = datetime.now()
        self.elapsed_times[taskstr] = str(self.finish_times[taskstr] - self.start_times[taskstr])
        
    def print_all_times(self):
        try:
            for task in self.tasks:
                if task not in self.elapsed_times:
                    self.end(task)
                if self.indents[task] > 0:
                    print("".join(["|", "---"*self.indents[task], ">"]),self.elapsed_times[task], task)
                else:
                    print(self.elapsed_times[task], task)
        except:
            print(self.elapsed_times)