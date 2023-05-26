import socket
import multiprocessing, os
import time
from datetime import datetime

End=b'\0'
def recv_end(conn):
    total_data=[]; data=''
    while True:
        data=conn.recv(1024)
        if End in data:
            total_data.append(data[:data.find(End)])
            break
        total_data.append(data)
    return b''.join(total_data)

def handleRequest(conn, addr):
   data = recv_end(conn)
   print('Process', os.getpid(), 'handling request from', addr, ':', data)
   time.sleep(10)

   if data == b'time':
       now = datetime.now()
       current_time = now.strftime("%H:%M:%S")
       conn.sendall(current_time.encode('utf-8') + b'\0')
   elif data == b'nuqneH?':
       conn.sendall(b'Qapla\'!\0')
   else:
       conn.sendall(data + b'\0')
   conn.close()
   
s = socket.socket()
host = socket.gethostname()
port = 12345
s.bind((host, port))

s.listen(5)
while True:
   conn, addr = s.accept()
   print('Process', os.getpid(),'got connection from', addr)
   process = multiprocessing.Process(target=handleRequest, args=(conn, addr))
   process.daemon = True
   process.start()

