from absl import app
from absl import flags

import zmq

import Utilities.GFP_Tools.nn_request_pb2 as nn_request

FLAGS = flags.FLAGS

flags.DEFINE_string("port", "", "Port tcp://host:port")
flags.DEFINE_boolean("shutdown", False, "Shut down the server")

def gfp_client(argv):
  context = zmq.Context()
  socket = context.socket(zmq.REQ)
  socket.connect(FLAGS.port)

  if FLAGS.shutdown:
    req = nn_request.Request()
    req.server_request.request = nn_request.ServerRequest.SHUTDOWN
    serialised = req.SerializeToString()
    socket.send(serialised)
    return

  for i in range(1):
    req = nn_request.Request()
    req.nn_request.smiles = "C1=NC(=CC2=C1C=CC=C2)NC(=O)NCC CHEMBL4076111"
    req.nn_request.id = "CHEMBL4076111"
    req.nn_request.nbrs = 10
    serialised = req.SerializeToString()
    socket.send(serialised)

    message = socket.recv()
    proto = nn_request.Reply()
    proto.ParseFromString(message);
    print(f"Received {proto}")


if __name__ == "__main__":
  app.run(gfp_client)
