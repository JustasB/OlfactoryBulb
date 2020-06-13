import os
from multiprocessing import Process
from time import sleep


class Blender:
    def __init__(self, keep=False, cmd_args=None, sleep=1):
        """
        Starts blender is a separate process. Use as `with Blender():`. Stops the process after end of `with`
        :param keep: When true, Blender process is not closed
        :param cmd_args: Any arguments to pass to the blender executable see: `blender --help` for options
        :param sleep: Number of seconds to delay further execution to allow Blender to finish initializing
        """
        self.keep = keep
        self.cmd_args = " " + ("" if cmd_args is None else cmd_args)
        self.sleep = sleep
        self.blender = Process(target=self.start_Blender)

    def start_Blender(self):
        os.system("blender"+self.cmd_args)

    def __enter__(self):
        self.blender.start()
        sleep(self.sleep)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.keep:
            os.system("pkill -9 blender") # This is a crude but potent way of stopping Blender
            self.blender.join()