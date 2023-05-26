import socket
import sys

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

message = b''

n = len(sys.argv)
if (n == 1):
    print("No message to server provided on the command line; ending program\n")
    exit(0)
else:
    for i in range(1, n):
        message += sys.argv[i].encode('UTF-8')
        if (i != n-1):
            message += b' '
            
message += End

print ("Sending message to server: " + message.decode('UTF-8') + "\n")

s = socket.socket()
host = socket.gethostname()
port = 12345

s.connect((host, port))

s.sendall(message)

print(recv_end(s))
s.close()
