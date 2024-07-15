

set -ex



pip check
py.test --pyargs zmq.tests.test_socket
exit 0
