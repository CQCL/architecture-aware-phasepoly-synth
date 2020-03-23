// Initial wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
// Resulting wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
cx q[7], q[3];
cx q[8], q[3];
cx q[9], q[3];
cx q[11], q[3];
cx q[13], q[3];
cx q[14], q[3];
cx q[15], q[3];
rz(0.25*pi) q[3];
cx q[7], q[3];
cx q[8], q[3];
cx q[9], q[3];
cx q[11], q[3];
cx q[13], q[3];
cx q[14], q[3];
cx q[15], q[3];
