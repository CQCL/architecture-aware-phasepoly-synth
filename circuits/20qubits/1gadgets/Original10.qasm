// Initial wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
// Resulting wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
cx q[2], q[0];
cx q[3], q[0];
cx q[6], q[0];
cx q[7], q[0];
cx q[9], q[0];
cx q[12], q[0];
cx q[13], q[0];
cx q[15], q[0];
cx q[17], q[0];
rz(0.25*pi) q[0];
cx q[2], q[0];
cx q[3], q[0];
cx q[6], q[0];
cx q[7], q[0];
cx q[9], q[0];
cx q[12], q[0];
cx q[13], q[0];
cx q[15], q[0];
cx q[17], q[0];
