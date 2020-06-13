import os, inspect
debug = False

class RunInClassDirectory:
    '''
    A utility class that temporarily switches the os working dir to the directory of the file in which the passed in
    class is located. The working dir is restored after the class is closed. The class should be used using the "with"
    statement e.g.:

    with RunInClassDirectory(NeuronCellClass):
        from neuron import h # <--- this will be executed in NeuronCellClass's folder
        etc...

    '''
    def __init__(self, the_class):
        self.the_class = the_class

    def __enter__(self):
        self.cur_dir = os.getcwd()
        class_file = inspect.getfile(self.the_class)
        class_dir = os.path.join(os.getcwd(), os.path.dirname(class_file))

        if debug:
            print('Temporarily changing directory to', class_dir)

        os.chdir(class_dir)

    def __exit__(self, exc_type, exc_value, traceback):
        if debug:
            print('Restoring working dir to', self.cur_dir)

        os.chdir(self.cur_dir)


class IsolatedCell(object):
    def close_window(self, name_contains):
        """
        Closes a NEURON window that matches the parameter string

        :param name_contains: The string to match in window name
        :return: Nothing
        """
        if not hasattr(self,"pwm"):
            from neuron import h
            self.pwm = h.PWManager()

        target_window_index = next(i
                                   for i in range(int(self.pwm.count()))
                                   if name_contains in self.pwm.name(i))

        self.pwm.close(target_window_index)
