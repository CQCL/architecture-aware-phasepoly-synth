// Initial wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
// Resulting wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
cx q[4], q[17];
cx q[2], q[17];
cx q[3], q[17];
cx q[5], q[17];
cx q[6], q[17];
cx q[8], q[17];
cx q[9], q[17];
cx q[10], q[17];
cx q[14], q[17];
cx q[18], q[17];
cx q[19], q[17];
rz(0.25*pi) q[17];
cx q[2], q[17];
cx q[0], q[17];
cx q[1], q[17];
cx q[3], q[17];
cx q[5], q[17];
cx q[7], q[17];
cx q[8], q[17];
cx q[10], q[17];
cx q[11], q[17];
cx q[14], q[17];
cx q[16], q[17];
rz(0.25*pi) q[17];
cx q[3], q[17];
cx q[18], q[17];
cx q[19], q[17];
cx q[5], q[17];
cx q[7], q[17];
cx q[8], q[17];
cx q[9], q[17];
cx q[10], q[17];
cx q[11], q[17];
cx q[12], q[17];
cx q[14], q[17];
cx q[16], q[17];
rz(0.25*pi) q[17];
cx q[5], q[17];
cx q[6], q[17];
cx q[7], q[17];
cx q[8], q[17];
cx q[9], q[17];
cx q[10], q[17];
cx q[11], q[17];
cx q[12], q[17];
cx q[13], q[17];
cx q[14], q[17];
cx q[15], q[17];
cx q[16], q[17];
rz(0.25*pi) q[17];
cx q[0], q[17];
cx q[1], q[17];
cx q[5], q[17];
cx q[6], q[17];
cx q[8], q[17];
cx q[9], q[17];
cx q[10], q[17];
cx q[11], q[17];
cx q[14], q[17];
cx q[15], q[17];
cx q[18], q[17];
cx q[19], q[17];
rz(0.25*pi) q[17];
cx q[4], q[17];
cx q[5], q[17];
cx q[13], q[17];
cx q[0], q[17];
cx q[3], q[17];
cx q[7], q[17];
cx q[10], q[17];
cx q[11], q[17];
cx q[14], q[17];
cx q[16], q[17];
cx q[18], q[17];
rz(0.25*pi) q[17];
cx q[0], q[17];
cx q[1], q[17];
cx q[18], q[17];
cx q[2], q[17];
cx q[6], q[17];
cx q[10], q[17];
cx q[14], q[17];
cx q[15], q[17];
cx q[16], q[17];
rz(0.25*pi) q[17];
cx q[2], q[17];
cx q[3], q[17];
cx q[6], q[17];
cx q[7], q[17];
cx q[8], q[17];
cx q[9], q[17];
cx q[10], q[17];
cx q[11], q[17];
cx q[14], q[17];
cx q[15], q[17];
cx q[16], q[17];
cx q[19], q[17];
rz(0.25*pi) q[17];
cx q[13], q[17];
cx q[1], q[17];
cx q[2], q[17];
cx q[3], q[17];
cx q[6], q[17];
cx q[7], q[17];
cx q[8], q[17];
cx q[10], q[17];
cx q[14], q[17];
cx q[16], q[17];
cx q[18], q[17];
cx q[19], q[17];
rz(0.25*pi) q[17];
cx q[19], q[17];
cx q[2], q[17];
cx q[3], q[17];
cx q[6], q[17];
cx q[8], q[17];
cx q[0], q[17];
cx q[10], q[17];
cx q[14], q[17];
cx q[16], q[17];
rz(0.25*pi) q[17];
cx q[0], q[17];
cx q[1], q[17];
cx q[7], q[17];
cx q[9], q[17];
cx q[10], q[17];
cx q[11], q[17];
cx q[12], q[17];
cx q[14], q[17];
cx q[15], q[17];
cx q[16], q[17];
cx q[18], q[17];
rz(0.25*pi) q[17];
cx q[2], q[17];
cx q[0], q[17];
cx q[1], q[17];
cx q[3], q[17];
cx q[7], q[17];
cx q[8], q[17];
cx q[9], q[17];
cx q[18], q[17];
rz(0.25*pi) q[17];
cx q[7], q[17];
cx q[9], q[17];
cx q[15], q[17];
cx q[0], q[17];
cx q[3], q[17];
cx q[6], q[17];
cx q[14], q[17];
cx q[1], q[17];
cx q[8], q[17];
cx q[12], q[17];
cx q[16], q[17];
rz(0.25*pi) q[17];
cx q[1], q[17];
cx q[8], q[17];
cx q[10], q[17];
cx q[11], q[17];
cx q[12], q[17];
cx q[16], q[17];
cx q[18], q[17];
rz(0.25*pi) q[17];
cx q[0], q[17];
cx q[1], q[17];
cx q[3], q[17];
cx q[6], q[17];
cx q[14], q[17];
cx q[16], q[17];
cx q[18], q[17];
rz(0.25*pi) q[17];
cx q[1], q[17];
cx q[5], q[17];
cx q[7], q[17];
cx q[12], q[17];
cx q[0], q[17];
cx q[6], q[17];
cx q[9], q[17];
cx q[16], q[17];
cx q[19], q[17];
rz(0.25*pi) q[17];
cx q[0], q[17];
cx q[6], q[17];
cx q[10], q[17];
cx q[15], q[17];
cx q[16], q[17];
cx q[19], q[17];
cx q[2], q[17];
cx q[3], q[17];
cx q[14], q[17];
rz(0.25*pi) q[17];
cx q[2], q[17];
cx q[3], q[17];
cx q[8], q[17];
cx q[9], q[17];
cx q[13], q[17];
cx q[14], q[17];
rz(0.25*pi) q[17];
cx q[5], q[8];
cx q[0], q[8];
cx q[2], q[8];
cx q[6], q[8];
cx q[7], q[8];
cx q[9], q[8];
cx q[11], q[8];
cx q[12], q[8];
cx q[14], q[8];
cx q[19], q[8];
rz(0.25*pi) q[8];
cx q[19], q[8];
cx q[9], q[8];
cx q[0], q[8];
cx q[1], q[8];
cx q[2], q[8];
cx q[7], q[8];
cx q[10], q[8];
cx q[12], q[8];
cx q[14], q[8];
cx q[16], q[8];
rz(0.25*pi) q[8];
cx q[1], q[8];
cx q[3], q[8];
cx q[4], q[8];
cx q[13], q[8];
cx q[14], q[8];
cx q[15], q[8];
cx q[2], q[8];
cx q[6], q[8];
cx q[7], q[8];
cx q[12], q[8];
cx q[16], q[8];
cx q[18], q[8];
rz(0.25*pi) q[8];
cx q[2], q[8];
cx q[6], q[8];
cx q[7], q[8];
cx q[10], q[8];
cx q[11], q[8];
cx q[12], q[8];
cx q[16], q[8];
cx q[18], q[8];
rz(0.25*pi) q[8];
cx q[0], q[8];
cx q[1], q[8];
cx q[3], q[8];
cx q[4], q[8];
cx q[6], q[8];
cx q[11], q[8];
cx q[12], q[8];
cx q[13], q[8];
cx q[14], q[8];
cx q[16], q[8];
cx q[18], q[8];
rz(0.25*pi) q[8];
cx q[2], q[8];
cx q[3], q[8];
cx q[6], q[8];
cx q[7], q[8];
cx q[9], q[8];
cx q[10], q[8];
cx q[11], q[8];
cx q[13], q[8];
cx q[14], q[8];
cx q[16], q[8];
rz(0.25*pi) q[8];
cx q[3], q[8];
cx q[5], q[8];
cx q[14], q[8];
cx q[0], q[8];
cx q[2], q[8];
cx q[11], q[8];
cx q[16], q[8];
cx q[19], q[8];
rz(0.25*pi) q[8];
cx q[0], q[8];
cx q[2], q[8];
cx q[4], q[8];
cx q[6], q[8];
cx q[7], q[8];
cx q[10], q[8];
cx q[11], q[8];
cx q[12], q[8];
cx q[13], q[8];
cx q[16], q[8];
cx q[19], q[8];
rz(0.25*pi) q[8];
cx q[0], q[16];
cx q[2], q[16];
cx q[3], q[16];
cx q[5], q[16];
cx q[6], q[16];
cx q[9], q[16];
cx q[18], q[16];
rz(0.25*pi) q[16];
cx q[2], q[16];
cx q[1], q[16];
cx q[4], q[16];
cx q[7], q[16];
cx q[9], q[16];
cx q[14], q[16];
cx q[18], q[16];
rz(0.25*pi) q[16];
cx q[9], q[16];
cx q[12], q[16];
cx q[3], q[16];
cx q[13], q[16];
cx q[5], q[16];
cx q[10], q[16];
cx q[14], q[16];
cx q[18], q[16];
cx q[19], q[16];
rz(0.25*pi) q[16];
cx q[5], q[16];
cx q[6], q[16];
cx q[10], q[16];
cx q[14], q[16];
cx q[15], q[16];
cx q[18], q[16];
cx q[19], q[16];
rz(0.25*pi) q[16];
cx q[3], q[16];
cx q[5], q[16];
cx q[10], q[16];
cx q[13], q[16];
cx q[18], q[16];
cx q[19], q[16];
rz(0.25*pi) q[16];
cx q[1], q[16];
cx q[4], q[16];
cx q[5], q[16];
cx q[6], q[16];
cx q[7], q[16];
cx q[10], q[16];
cx q[11], q[16];
cx q[13], q[16];
cx q[15], q[16];
cx q[19], q[16];
rz(0.25*pi) q[16];
cx q[0], q[16];
cx q[5], q[16];
cx q[9], q[16];
cx q[14], q[16];
cx q[15], q[16];
cx q[18], q[16];
cx q[19], q[16];
cx q[1], q[16];
cx q[2], q[16];
cx q[3], q[16];
cx q[4], q[16];
cx q[12], q[16];
rz(0.25*pi) q[16];
cx q[1], q[16];
cx q[2], q[16];
cx q[3], q[16];
cx q[4], q[16];
cx q[6], q[16];
cx q[10], q[16];
cx q[12], q[16];
cx q[13], q[16];
rz(0.25*pi) q[16];
cx q[5], q[14];
cx q[9], q[14];
cx q[0], q[14];
cx q[2], q[14];
cx q[13], q[14];
cx q[19], q[14];
cx q[1], q[14];
cx q[4], q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[12], q[14];
cx q[15], q[14];
cx q[18], q[14];
rz(0.25*pi) q[14];
cx q[1], q[14];
cx q[4], q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[10], q[14];
cx q[12], q[14];
cx q[15], q[14];
cx q[18], q[14];
rz(0.25*pi) q[14];
cx q[0], q[14];
cx q[2], q[14];
cx q[4], q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[10], q[14];
cx q[11], q[14];
cx q[12], q[14];
cx q[13], q[14];
cx q[19], q[14];
rz(0.25*pi) q[14];
cx q[9], q[14];
cx q[3], q[14];
cx q[5], q[14];
cx q[0], q[14];
cx q[1], q[14];
cx q[2], q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[10], q[14];
cx q[12], q[14];
cx q[19], q[14];
rz(0.25*pi) q[14];
cx q[2], q[14];
cx q[4], q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[10], q[14];
cx q[12], q[14];
cx q[13], q[14];
cx q[15], q[14];
cx q[18], q[14];
cx q[19], q[14];
rz(0.25*pi) q[14];
cx q[1], q[14];
cx q[0], q[14];
cx q[2], q[14];
cx q[12], q[14];
cx q[18], q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[10], q[14];
cx q[11], q[14];
cx q[15], q[14];
rz(0.25*pi) q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[10], q[14];
cx q[11], q[14];
cx q[15], q[14];
cx q[19], q[14];
rz(0.25*pi) q[14];
cx q[2], q[14];
cx q[4], q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[10], q[14];
cx q[12], q[14];
cx q[13], q[14];
cx q[15], q[14];
cx q[18], q[14];
rz(0.25*pi) q[14];
cx q[0], q[14];
cx q[2], q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[11], q[14];
cx q[12], q[14];
cx q[15], q[14];
cx q[19], q[14];
rz(0.25*pi) q[14];
cx q[2], q[14];
cx q[4], q[14];
cx q[5], q[14];
cx q[6], q[14];
cx q[7], q[14];
cx q[15], q[14];
cx q[19], q[14];
rz(0.25*pi) q[14];
cx q[1], q[14];
cx q[3], q[14];
cx q[4], q[14];
cx q[5], q[14];
cx q[10], q[14];
cx q[12], q[14];
cx q[13], q[14];
cx q[19], q[14];
rz(0.25*pi) q[14];
cx q[5], q[1];
cx q[0], q[1];
cx q[2], q[1];
cx q[3], q[1];
cx q[9], q[1];
cx q[10], q[1];
cx q[11], q[1];
cx q[18], q[1];
cx q[19], q[1];
rz(0.25*pi) q[1];
cx q[0], q[1];
cx q[11], q[1];
cx q[2], q[1];
cx q[3], q[1];
cx q[12], q[1];
cx q[13], q[1];
cx q[19], q[1];
cx q[4], q[1];
cx q[7], q[1];
cx q[9], q[1];
cx q[18], q[1];
rz(0.25*pi) q[1];
cx q[4], q[1];
cx q[7], q[1];
cx q[9], q[1];
cx q[10], q[1];
cx q[18], q[1];
rz(0.25*pi) q[1];
cx q[2], q[1];
cx q[3], q[1];
cx q[4], q[1];
cx q[6], q[1];
cx q[7], q[1];
cx q[9], q[1];
cx q[12], q[1];
cx q[13], q[1];
cx q[19], q[1];
rz(0.25*pi) q[1];
cx q[3], q[0];
cx q[4], q[0];
cx q[6], q[0];
cx q[7], q[0];
cx q[9], q[0];
cx q[15], q[0];
cx q[19], q[0];
rz(0.25*pi) q[0];
cx q[2], q[1];
cx q[3], q[1];
cx q[3], q[0];
cx q[4], q[1];
cx q[4], q[0];
cx q[5], q[1];
cx q[6], q[1];
cx q[6], q[0];
cx q[7], q[1];
cx q[7], q[0];
cx q[9], q[8];
cx q[9], q[0];
cx q[10], q[8];
cx q[11], q[8];
cx q[12], q[8];
cx q[16], q[8];
cx q[15], q[0];
cx q[8], q[1];
cx q[18], q[8];
cx q[19], q[8];
cx q[19], q[0];
cx q[8], q[1];
cx q[18], q[14];
cx q[19], q[16];
cx q[15], q[17];
cx q[15], q[16];
cx q[14], q[17];
cx q[13], q[17];
cx q[13], q[14];
cx q[12], q[17];
cx q[12], q[16];
cx q[12], q[14];
cx q[11], q[16];
cx q[10], q[17];
cx q[10], q[16];
cx q[6], q[17];
cx q[5], q[17];
cx q[4], q[17];
cx q[4], q[14];
cx q[4], q[8];
cx q[3], q[8];
cx q[2], q[17];
cx q[2], q[8];
cx q[1], q[17];
cx q[1], q[14];
cx q[1], q[8];
cx q[0], q[14];
cx q[0], q[8];
cx q[3], q[16];
