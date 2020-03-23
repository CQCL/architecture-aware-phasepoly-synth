// Initial wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
// Resulting wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[5], q[4];
cx q[6], q[4];
cx q[7], q[4];
cx q[8], q[4];
cx q[13], q[4];
cx q[15], q[4];
rz(0.25*pi) q[4];
cx q[5], q[4];
cx q[6], q[4];
cx q[7], q[4];
cx q[8], q[4];
cx q[13], q[4];
cx q[15], q[4];
