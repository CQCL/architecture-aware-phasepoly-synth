// Initial wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
// Resulting wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[1], q[0];
cx q[2], q[0];
cx q[5], q[0];
cx q[8], q[0];
cx q[11], q[0];
cx q[13], q[0];
rz(0.25*pi) q[0];
cx q[1], q[0];
cx q[2], q[0];
cx q[5], q[0];
cx q[6], q[0];
cx q[8], q[0];
cx q[9], q[0];
cx q[11], q[0];
cx q[13], q[0];
rz(0.25*pi) q[0];
cx q[2], q[9];
cx q[4], q[9];
cx q[5], q[9];
cx q[7], q[9];
cx q[8], q[9];
cx q[11], q[9];
cx q[12], q[9];
cx q[15], q[9];
cx q[6], q[9];
rz(0.25*pi) q[9];
cx q[6], q[9];
cx q[10], q[9];
rz(0.25*pi) q[9];
cx q[1], q[2];
cx q[3], q[2];
cx q[4], q[2];
cx q[7], q[2];
cx q[8], q[2];
cx q[10], q[2];
cx q[11], q[2];
cx q[12], q[2];
cx q[13], q[2];
rz(0.25*pi) q[2];
cx q[1], q[3];
cx q[4], q[3];
cx q[6], q[3];
cx q[7], q[3];
cx q[11], q[3];
cx q[14], q[3];
rz(0.25*pi) q[3];
cx q[7], q[3];
cx q[11], q[3];
cx q[1], q[3];
cx q[6], q[3];
cx q[4], q[3];
cx q[5], q[3];
cx q[15], q[3];
rz(0.25*pi) q[3];
cx q[4], q[3];
cx q[5], q[3];
cx q[8], q[3];
cx q[10], q[3];
cx q[13], q[3];
cx q[14], q[3];
cx q[15], q[3];
rz(0.25*pi) q[3];
cx q[1], q[3];
cx q[4], q[3];
cx q[5], q[3];
cx q[6], q[3];
cx q[8], q[3];
cx q[13], q[3];
cx q[14], q[3];
cx q[15], q[3];
rz(0.25*pi) q[3];
cx q[5], q[1];
cx q[7], q[1];
cx q[11], q[1];
cx q[13], q[1];
cx q[14], q[1];
rz(0.25*pi) q[1];
cx q[1], q[0];
cx q[2], q[0];
cx q[3], q[2];
cx q[3], q[0];
cx q[4], q[2];
cx q[5], q[3];
cx q[5], q[2];
cx q[5], q[1];
cx q[5], q[0];
cx q[6], q[3];
cx q[6], q[2];
cx q[7], q[2];
cx q[7], q[1];
cx q[7], q[0];
cx q[8], q[2];
cx q[9], q[0];
cx q[10], q[9];
cx q[10], q[3];
cx q[10], q[0];
cx q[11], q[9];
cx q[11], q[2];
cx q[11], q[1];
cx q[11], q[0];
cx q[12], q[9];
cx q[12], q[2];
cx q[13], q[2];
cx q[13], q[1];
cx q[14], q[3];
cx q[14], q[2];
cx q[14], q[1];
cx q[15], q[9];
cx q[15], q[3];
cx q[15], q[2];
cx q[1], q[0];
cx q[8], q[9];
cx q[7], q[9];
cx q[5], q[9];
cx q[4], q[9];
cx q[2], q[9];
cx q[1], q[3];
