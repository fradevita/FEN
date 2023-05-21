import sys
import os

SOURCE=sys.argv[1]

command = "grep 'use' " + SOURCE + " | awk '{print $2}'"
stream = os.popen(command)
list = stream.read()
object = ''
object_list = [""]
counter = 0
for c in list:
    if (c == '\n'):
        if (counter == 0):
            object_list[0] = object
        else:
            object_list.append(object)
        object = ''
        counter = counter + 1
    else:
        object = object + c

# sostituisco '_mod' con '.o'
f = open('.objects', 'w')
for obj in object_list:
    if obj != 'mpi':
        f.write(obj.replace("_mod", ".o")+'\n')
        print(obj.replace("_mod", ".o")+'\n')

f.close()