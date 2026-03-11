pyenv('Version','C:\Users\jwooj\AppData\Local\Programs\Python\Python312\python.exe');

pyModule = py.importlib.import_module('numpy');

pyArray = pyModule.load('pointcloud_frame.npy');
matlabArray = double(pyArray);