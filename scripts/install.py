import platform
import os
import os.path
import shutil
import subprocess

scripts_path = os.path.realpath(os.path.dirname(__file__))

root_path = os.path.abspath(f"{scripts_path}\..")

src_path = os.path.abspath(f"{root_path}\src")
obj_path = os.path.abspath(f"{root_path}\obj")
libs_path = os.path.abspath(f"{root_path}\libs")
dll_path = os.path.abspath(f"{root_path}\DLL")

def run_command(cmd):
  p = subprocess.Popen(cmd);  
  res = p.communicate()
  
  retcode = p.returncode
  if retcode != 0:
    raise subprocess.CalledProcessError(returncode=retcode, cmd=" ".join(cmd), stderr=res[1])
    
def windows_create_obj_file(source_name, obj_name):
    src_file = os.path.join(src_path, f"{source_name}.c")
    obj_file = os.path.join(obj_path, f"{obj_name}.o")
    cmd = ["gcc", "-c", "-o", obj_file, src_file]
    run_command(cmd)
    


def create_dependencies_args(dependencies):
  args = map(lambda x: "-l{}".format(x), dependencies)
  return " ".join(args)
    
def windows_create_lib_file(obj_name, dependencies):
    obj_file = os.path.join(obj_path, f"{obj_name}.o")
    dll_file = os.path.join(dll_path, f"{obj_name}.dll")
    libs_file = os.path.join(libs_path, f"lib{obj_name}.a")
    if len(dependencies) == 0:
      cmd = "gcc -o {} {} -s -shared -Wl,--subsystem,windows,--out-implib,{}".format(
        dll_file, obj_file, libs_file)
      
      run_command(cmd)
    
    else:
      dependencies_args = create_dependencies_args(dependencies)
      cmd = "gcc -o {} {} -s -shared -Wl,--subsystem,windows,--out-implib,{} -L{} {}".format(
        dll_file, obj_file, libs_file, libs_path, dependencies_args)

      run_command(cmd)

if __name__ == "__main__":
    if platform.system() == "Windows":
        os.mkdir(obj_path)
        os.mkdir(libs_path)

        if not os.path.isdir(dll_path):
            os.mkdir(dll_path)

        source = [("complex", "complex"), ("fft", "fft"),
                  ("iterative_fft", "itfft"), ("conv", "conv")]

        dependencies = [[], ["complex"], ["fft", "complex"], ["itfft", "fft", "complex"]]
        
        try:
          for (source_name, obj_name), dependency in zip(source, dependencies):   
            windows_create_obj_file(source_name, obj_name)
            windows_create_lib_file(obj_name, dependency)
            
        except Exception as e:
          print(e)
          
        finally:
          shutil.rmtree(obj_path)
          shutil.rmtree(libs_path)
          pass
        
        print("OK")

    else:
        os.system("./scripts.sh")
