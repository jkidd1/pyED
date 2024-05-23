import os
from contextlib import contextmanager

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

class reader:
        def __init__(self, name, strip_lines = True):
                self.filename = name
                self.exists = self.isFile()
                if self.exists:
                        self.contents = self.get_lines(strip_lines)
                else:
                        self.contents = None
        def isFile(self):
                from os import listdir
                from os.path import isfile, join
                directory = '.'
                if '/' in self.filename:
                        directory = '/'.join(self.filename.split('/')[:-1])
                dir_files = [f for f in listdir(directory) if isfile(join(directory, f))]
                if self.filename.split('/')[-1] in dir_files:
                        return True
                return False
        def get_lines(self, strip_lines):
            with open(self.filename, 'r') as file:
                if strip_lines:
                    return [line.strip() for line in file.readlines()]
                return [line.replace('\n', '') for line in file.readlines()]

class writer:
        def __init__(self, name, strip = False):
                self.filename = name
                self.strip = strip
                self.exists = reader(self.filename).exists
        def replace(self, old_string, new_string):
            old_file = reader(self.filename, self.strip)
            if old_file:
                    old_lines = old_file.contents
                    new_lines = '\n'.join([line.replace(old_string, new_string) for line in old_lines])
                    with open(self.filename, 'w') as new_file:
                        new_file.write(new_lines)
