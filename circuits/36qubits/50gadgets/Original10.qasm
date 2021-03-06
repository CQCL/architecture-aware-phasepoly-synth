// Initial wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
// Resulting wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
OPENQASM 2.0;
include "qelib1.inc";
qreg q[36];
cx q[0], q[15];
cx q[13], q[15];
cx q[5], q[15];
cx q[2], q[15];
cx q[22], q[15];
cx q[29], q[15];
cx q[3], q[15];
cx q[4], q[15];
cx q[7], q[15];
cx q[9], q[15];
cx q[14], q[15];
cx q[17], q[15];
cx q[20], q[15];
cx q[24], q[15];
cx q[26], q[15];
cx q[27], q[15];
cx q[32], q[15];
cx q[34], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[3], q[15];
cx q[4], q[15];
cx q[14], q[15];
cx q[26], q[15];
cx q[1], q[15];
cx q[6], q[15];
cx q[16], q[15];
cx q[25], q[15];
cx q[32], q[15];
cx q[34], q[15];
cx q[7], q[15];
cx q[8], q[15];
cx q[12], q[15];
cx q[20], q[15];
cx q[28], q[15];
cx q[31], q[15];
rz(0.25*pi) q[15];
cx q[7], q[15];
cx q[8], q[15];
cx q[9], q[15];
cx q[11], q[15];
cx q[12], q[15];
cx q[17], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[23], q[15];
cx q[24], q[15];
cx q[27], q[15];
cx q[28], q[15];
cx q[31], q[15];
cx q[33], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[1], q[15];
cx q[6], q[15];
cx q[8], q[15];
cx q[11], q[15];
cx q[12], q[15];
cx q[16], q[15];
cx q[17], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[25], q[15];
cx q[27], q[15];
cx q[32], q[15];
cx q[34], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[22], q[15];
cx q[11], q[15];
cx q[10], q[15];
cx q[1], q[15];
cx q[3], q[15];
cx q[6], q[15];
cx q[12], q[15];
cx q[16], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[24], q[15];
cx q[26], q[15];
cx q[27], q[15];
cx q[28], q[15];
cx q[29], q[15];
cx q[31], q[15];
cx q[33], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[24], q[15];
cx q[4], q[15];
cx q[2], q[15];
cx q[8], q[15];
cx q[9], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[17], q[15];
cx q[18], q[15];
cx q[20], q[15];
cx q[25], q[15];
cx q[26], q[15];
cx q[29], q[15];
cx q[30], q[15];
cx q[32], q[15];
cx q[33], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[9], q[15];
cx q[23], q[15];
cx q[33], q[15];
cx q[34], q[15];
cx q[35], q[15];
cx q[1], q[15];
cx q[12], q[15];
cx q[25], q[15];
cx q[2], q[15];
cx q[6], q[15];
cx q[7], q[15];
cx q[16], q[15];
cx q[19], q[15];
cx q[21], q[15];
cx q[27], q[15];
cx q[29], q[15];
cx q[31], q[15];
cx q[32], q[15];
rz(0.25*pi) q[15];
cx q[6], q[15];
cx q[7], q[15];
cx q[8], q[15];
cx q[16], q[15];
cx q[17], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[21], q[15];
cx q[27], q[15];
cx q[29], q[15];
cx q[30], q[15];
cx q[31], q[15];
cx q[32], q[15];
rz(0.25*pi) q[15];
cx q[2], q[15];
cx q[3], q[15];
cx q[7], q[15];
cx q[8], q[15];
cx q[14], q[15];
cx q[16], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[26], q[15];
cx q[27], q[15];
cx q[28], q[15];
cx q[30], q[15];
cx q[32], q[15];
rz(0.25*pi) q[15];
cx q[1], q[15];
cx q[6], q[15];
cx q[7], q[15];
cx q[8], q[15];
cx q[12], q[15];
cx q[17], q[15];
cx q[18], q[15];
cx q[20], q[15];
cx q[25], q[15];
cx q[26], q[15];
cx q[30], q[15];
rz(0.25*pi) q[15];
cx q[1], q[15];
cx q[3], q[15];
cx q[4], q[15];
cx q[9], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[17], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[28], q[15];
cx q[34], q[15];
rz(0.25*pi) q[15];
cx q[1], q[15];
cx q[2], q[15];
cx q[4], q[15];
cx q[6], q[15];
cx q[7], q[15];
cx q[9], q[15];
cx q[10], q[15];
cx q[16], q[15];
cx q[17], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[23], q[15];
cx q[25], q[15];
cx q[26], q[15];
cx q[27], q[15];
cx q[28], q[15];
cx q[29], q[15];
cx q[30], q[15];
cx q[33], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[1], q[15];
cx q[3], q[15];
cx q[8], q[15];
cx q[11], q[15];
cx q[16], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[27], q[15];
cx q[31], q[15];
cx q[2], q[15];
cx q[7], q[15];
cx q[12], q[15];
cx q[23], q[15];
cx q[24], q[15];
cx q[25], q[15];
cx q[32], q[15];
cx q[33], q[15];
rz(0.25*pi) q[15];
cx q[2], q[15];
cx q[4], q[15];
cx q[7], q[15];
cx q[9], q[15];
cx q[12], q[15];
cx q[19], q[15];
cx q[23], q[15];
cx q[24], q[15];
cx q[25], q[15];
cx q[28], q[15];
cx q[30], q[15];
cx q[32], q[15];
cx q[33], q[15];
cx q[34], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[4], q[15];
cx q[5], q[15];
cx q[16], q[15];
cx q[17], q[15];
cx q[1], q[15];
cx q[6], q[15];
cx q[7], q[15];
cx q[8], q[15];
cx q[10], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[21], q[15];
cx q[22], q[15];
cx q[23], q[15];
cx q[25], q[15];
cx q[27], q[15];
cx q[28], q[15];
cx q[34], q[15];
rz(0.25*pi) q[15];
cx q[6], q[15];
cx q[7], q[15];
cx q[8], q[15];
cx q[30], q[15];
cx q[34], q[15];
cx q[1], q[15];
cx q[2], q[15];
cx q[3], q[15];
cx q[24], q[15];
cx q[28], q[15];
cx q[29], q[15];
cx q[10], q[15];
cx q[11], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[22], q[15];
cx q[25], q[15];
cx q[26], q[15];
rz(0.25*pi) q[15];
cx q[10], q[15];
cx q[11], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[22], q[15];
cx q[23], q[15];
cx q[25], q[15];
cx q[26], q[15];
cx q[27], q[15];
cx q[31], q[15];
cx q[33], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[1], q[15];
cx q[2], q[15];
cx q[3], q[15];
cx q[9], q[15];
cx q[10], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[18], q[15];
cx q[20], q[15];
cx q[22], q[15];
cx q[24], q[15];
cx q[25], q[15];
cx q[27], q[15];
cx q[28], q[15];
cx q[29], q[15];
cx q[31], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[3], q[15];
cx q[6], q[15];
cx q[7], q[15];
cx q[8], q[15];
cx q[11], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[17], q[15];
cx q[18], q[15];
cx q[21], q[15];
cx q[24], q[15];
cx q[25], q[15];
cx q[28], q[15];
cx q[31], q[15];
cx q[34], q[15];
rz(0.25*pi) q[15];
cx q[13], q[15];
cx q[35], q[15];
cx q[1], q[15];
cx q[5], q[15];
cx q[7], q[15];
cx q[9], q[15];
cx q[11], q[15];
cx q[17], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[21], q[15];
cx q[24], q[15];
cx q[26], q[15];
cx q[27], q[15];
cx q[28], q[15];
cx q[30], q[15];
cx q[31], q[15];
cx q[32], q[15];
cx q[33], q[15];
rz(0.25*pi) q[15];
cx q[1], q[15];
cx q[24], q[15];
cx q[28], q[15];
cx q[5], q[15];
cx q[7], q[15];
cx q[25], q[15];
cx q[27], q[15];
cx q[2], q[15];
cx q[4], q[15];
cx q[10], q[15];
cx q[16], q[15];
cx q[17], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[22], q[15];
cx q[26], q[15];
cx q[30], q[15];
cx q[34], q[15];
rz(0.25*pi) q[15];
cx q[2], q[15];
cx q[6], q[15];
cx q[8], q[15];
cx q[9], q[15];
cx q[10], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[32], q[15];
cx q[3], q[15];
cx q[4], q[15];
cx q[14], q[15];
cx q[16], q[15];
cx q[17], q[15];
cx q[22], q[15];
cx q[23], q[15];
cx q[30], q[15];
cx q[33], q[15];
rz(0.25*pi) q[15];
cx q[3], q[15];
cx q[4], q[15];
cx q[11], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[16], q[15];
cx q[17], q[15];
cx q[22], q[15];
cx q[23], q[15];
cx q[26], q[15];
cx q[29], q[15];
cx q[30], q[15];
cx q[31], q[15];
cx q[33], q[15];
cx q[34], q[15];
rz(0.25*pi) q[15];
cx q[2], q[15];
cx q[3], q[15];
cx q[4], q[15];
cx q[5], q[15];
cx q[7], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[21], q[15];
cx q[22], q[15];
cx q[25], q[15];
cx q[27], q[15];
cx q[29], q[15];
cx q[30], q[15];
rz(0.25*pi) q[15];
cx q[0], q[15];
cx q[8], q[15];
cx q[2], q[15];
cx q[4], q[15];
cx q[6], q[15];
cx q[7], q[15];
cx q[10], q[15];
cx q[13], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[22], q[15];
cx q[23], q[15];
cx q[24], q[15];
cx q[25], q[15];
cx q[27], q[15];
cx q[30], q[15];
cx q[33], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[7], q[15];
cx q[1], q[15];
cx q[2], q[15];
cx q[4], q[15];
cx q[10], q[15];
cx q[12], q[15];
cx q[18], q[15];
cx q[28], q[15];
cx q[30], q[15];
cx q[31], q[15];
cx q[3], q[15];
cx q[5], q[15];
cx q[16], q[15];
cx q[25], q[15];
rz(0.25*pi) q[15];
cx q[3], q[15];
cx q[5], q[15];
cx q[9], q[15];
cx q[14], q[15];
cx q[16], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[22], q[15];
cx q[24], q[15];
cx q[25], q[15];
cx q[26], q[15];
cx q[29], q[15];
cx q[32], q[15];
cx q[33], q[15];
cx q[34], q[15];
rz(0.25*pi) q[15];
cx q[30], q[15];
cx q[2], q[15];
cx q[3], q[15];
cx q[4], q[15];
cx q[6], q[15];
cx q[9], q[15];
cx q[13], q[15];
cx q[16], q[15];
cx q[17], q[15];
cx q[19], q[15];
cx q[21], q[15];
cx q[23], q[15];
cx q[24], q[15];
cx q[25], q[15];
cx q[26], q[15];
cx q[28], q[15];
cx q[32], q[15];
cx q[33], q[15];
cx q[34], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[34], q[15];
cx q[1], q[15];
cx q[3], q[15];
cx q[18], q[15];
cx q[2], q[15];
cx q[5], q[15];
cx q[9], q[15];
cx q[12], q[15];
cx q[13], q[15];
cx q[14], q[15];
cx q[16], q[15];
cx q[20], q[15];
cx q[27], q[15];
cx q[28], q[15];
cx q[32], q[15];
cx q[33], q[15];
rz(0.25*pi) q[15];
cx q[5], q[15];
cx q[6], q[15];
cx q[9], q[15];
cx q[25], q[15];
cx q[31], q[15];
cx q[2], q[15];
cx q[4], q[15];
cx q[13], q[15];
cx q[17], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[29], q[15];
cx q[10], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[16], q[15];
cx q[23], q[15];
cx q[26], q[15];
cx q[32], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[10], q[15];
cx q[11], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[16], q[15];
cx q[21], q[15];
cx q[23], q[15];
cx q[24], q[15];
cx q[26], q[15];
cx q[28], q[15];
cx q[32], q[15];
cx q[33], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[2], q[15];
cx q[4], q[15];
cx q[10], q[15];
cx q[13], q[15];
cx q[16], q[15];
cx q[17], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[22], q[15];
cx q[26], q[15];
cx q[27], q[15];
cx q[29], q[15];
cx q[32], q[15];
cx q[35], q[15];
rz(0.25*pi) q[15];
cx q[2], q[15];
cx q[3], q[15];
cx q[6], q[15];
cx q[12], q[15];
cx q[13], q[15];
cx q[16], q[15];
cx q[18], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[23], q[15];
cx q[24], q[15];
cx q[26], q[15];
cx q[27], q[15];
cx q[28], q[15];
cx q[32], q[15];
rz(0.25*pi) q[15];
cx q[1], q[15];
cx q[2], q[15];
cx q[4], q[15];
cx q[5], q[15];
cx q[9], q[15];
cx q[10], q[15];
cx q[11], q[15];
cx q[12], q[15];
cx q[14], q[15];
cx q[21], q[15];
cx q[23], q[15];
cx q[25], q[15];
cx q[26], q[15];
rz(0.25*pi) q[15];
cx q[1], q[15];
cx q[4], q[15];
cx q[7], q[15];
cx q[8], q[15];
cx q[9], q[15];
cx q[10], q[15];
cx q[11], q[15];
cx q[16], q[15];
cx q[19], q[15];
cx q[20], q[15];
cx q[21], q[15];
cx q[22], q[15];
cx q[23], q[15];
cx q[25], q[15];
cx q[26], q[15];
cx q[27], q[15];
cx q[29], q[15];
cx q[30], q[15];
cx q[33], q[15];
cx q[34], q[15];
rz(0.25*pi) q[15];
cx q[6], q[30];
cx q[33], q[30];
cx q[0], q[30];
cx q[13], q[30];
cx q[23], q[30];
cx q[24], q[30];
cx q[25], q[30];
cx q[31], q[30];
cx q[1], q[30];
cx q[3], q[30];
cx q[18], q[30];
cx q[27], q[30];
rz(0.25*pi) q[30];
cx q[1], q[30];
cx q[2], q[30];
cx q[3], q[30];
cx q[5], q[30];
cx q[7], q[30];
cx q[8], q[30];
cx q[10], q[30];
cx q[11], q[30];
cx q[12], q[30];
cx q[14], q[30];
cx q[17], q[30];
cx q[18], q[30];
cx q[26], q[30];
cx q[27], q[30];
rz(0.25*pi) q[30];
cx q[0], q[30];
cx q[1], q[30];
cx q[4], q[30];
cx q[5], q[30];
cx q[7], q[30];
cx q[9], q[30];
cx q[11], q[30];
cx q[13], q[30];
cx q[14], q[30];
cx q[16], q[30];
cx q[19], q[30];
cx q[20], q[30];
cx q[22], q[30];
cx q[23], q[30];
cx q[24], q[30];
cx q[25], q[30];
cx q[26], q[30];
cx q[31], q[30];
rz(0.25*pi) q[30];
cx q[0], q[3];
cx q[1], q[3];
cx q[5], q[3];
cx q[6], q[3];
cx q[7], q[3];
cx q[10], q[3];
cx q[11], q[3];
cx q[13], q[3];
cx q[24], q[3];
cx q[27], q[3];
cx q[34], q[3];
cx q[2], q[3];
cx q[14], q[3];
cx q[17], q[3];
cx q[19], q[3];
cx q[21], q[3];
cx q[22], q[3];
cx q[25], q[3];
cx q[28], q[3];
cx q[33], q[3];
cx q[35], q[3];
rz(0.25*pi) q[3];
cx q[2], q[3];
cx q[12], q[3];
cx q[14], q[3];
cx q[16], q[3];
cx q[17], q[3];
cx q[19], q[3];
cx q[20], q[3];
cx q[21], q[3];
cx q[22], q[3];
cx q[23], q[3];
cx q[25], q[3];
cx q[26], q[3];
cx q[28], q[3];
cx q[29], q[3];
cx q[33], q[3];
cx q[35], q[3];
rz(0.25*pi) q[3];
cx q[1], q[3];
cx q[14], q[3];
cx q[0], q[3];
cx q[2], q[3];
cx q[4], q[3];
cx q[5], q[3];
cx q[11], q[3];
cx q[13], q[3];
cx q[17], q[3];
cx q[18], q[3];
cx q[19], q[3];
cx q[20], q[3];
cx q[21], q[3];
cx q[24], q[3];
cx q[25], q[3];
cx q[27], q[3];
cx q[28], q[3];
cx q[32], q[3];
cx q[33], q[3];
rz(0.25*pi) q[3];
cx q[20], q[3];
cx q[33], q[3];
cx q[7], q[3];
cx q[19], q[3];
cx q[0], q[3];
cx q[4], q[3];
cx q[2], q[3];
cx q[11], q[3];
cx q[12], q[3];
cx q[13], q[3];
cx q[5], q[3];
cx q[10], q[3];
cx q[16], q[3];
cx q[17], q[3];
cx q[18], q[3];
cx q[21], q[3];
cx q[22], q[3];
cx q[28], q[3];
rz(0.25*pi) q[3];
cx q[5], q[3];
cx q[9], q[3];
cx q[10], q[3];
cx q[16], q[3];
cx q[17], q[3];
cx q[18], q[3];
cx q[21], q[3];
cx q[22], q[3];
cx q[24], q[3];
cx q[25], q[3];
cx q[26], q[3];
cx q[28], q[3];
cx q[29], q[3];
cx q[34], q[3];
cx q[35], q[3];
rz(0.25*pi) q[3];
cx q[2], q[3];
cx q[6], q[3];
cx q[11], q[3];
cx q[12], q[3];
cx q[13], q[3];
cx q[18], q[3];
cx q[21], q[3];
cx q[22], q[3];
cx q[23], q[3];
cx q[24], q[3];
cx q[25], q[3];
cx q[26], q[3];
cx q[27], q[3];
cx q[29], q[3];
cx q[31], q[3];
cx q[35], q[3];
rz(0.25*pi) q[3];
cx q[0], q[3];
cx q[4], q[3];
cx q[8], q[3];
cx q[10], q[3];
cx q[11], q[3];
cx q[12], q[3];
cx q[16], q[3];
cx q[22], q[3];
cx q[23], q[3];
cx q[24], q[3];
cx q[27], q[3];
cx q[32], q[3];
rz(0.25*pi) q[3];
cx q[2], q[3];
cx q[5], q[3];
cx q[7], q[3];
cx q[16], q[3];
cx q[17], q[3];
cx q[18], q[3];
cx q[19], q[3];
cx q[21], q[3];
cx q[22], q[3];
cx q[24], q[3];
cx q[27], q[3];
cx q[28], q[3];
cx q[29], q[3];
cx q[31], q[3];
rz(0.25*pi) q[3];
cx q[5], q[3];
cx q[6], q[3];
cx q[7], q[3];
cx q[8], q[3];
cx q[9], q[3];
cx q[12], q[3];
cx q[14], q[3];
cx q[17], q[3];
cx q[18], q[3];
cx q[19], q[3];
cx q[21], q[3];
cx q[22], q[3];
cx q[25], q[3];
cx q[26], q[3];
cx q[32], q[3];
cx q[33], q[3];
cx q[34], q[3];
rz(0.25*pi) q[3];
cx q[8], q[0];
cx q[10], q[0];
cx q[12], q[0];
cx q[13], q[0];
cx q[14], q[0];
cx q[16], q[0];
cx q[17], q[0];
cx q[18], q[0];
cx q[19], q[0];
cx q[24], q[0];
cx q[29], q[0];
cx q[31], q[0];
cx q[32], q[0];
cx q[33], q[0];
rz(0.25*pi) q[0];
cx q[1], q[5];
cx q[2], q[5];
cx q[4], q[5];
cx q[6], q[5];
cx q[7], q[5];
cx q[9], q[5];
cx q[10], q[5];
cx q[12], q[5];
cx q[14], q[5];
cx q[16], q[5];
cx q[17], q[5];
cx q[20], q[5];
cx q[22], q[5];
cx q[24], q[5];
cx q[26], q[5];
cx q[28], q[5];
cx q[33], q[5];
cx q[34], q[5];
cx q[35], q[5];
rz(0.25*pi) q[5];
cx q[2], q[1];
cx q[4], q[1];
cx q[6], q[1];
cx q[8], q[1];
cx q[9], q[1];
cx q[11], q[1];
cx q[12], q[1];
cx q[13], q[1];
cx q[19], q[1];
cx q[22], q[1];
cx q[24], q[1];
cx q[28], q[1];
cx q[29], q[1];
cx q[34], q[1];
cx q[35], q[1];
rz(0.25*pi) q[1];
cx q[2], q[1];
cx q[4], q[3];
cx q[4], q[1];
cx q[6], q[5];
cx q[6], q[3];
cx q[6], q[1];
cx q[7], q[5];
cx q[8], q[1];
cx q[8], q[0];
cx q[9], q[5];
cx q[9], q[1];
cx q[10], q[5];
cx q[10], q[0];
cx q[11], q[3];
cx q[11], q[1];
cx q[12], q[5];
cx q[12], q[3];
cx q[12], q[1];
cx q[12], q[0];
cx q[13], q[1];
cx q[13], q[0];
cx q[14], q[5];
cx q[14], q[0];
cx q[16], q[15];
cx q[16], q[5];
cx q[16], q[3];
cx q[16], q[0];
cx q[17], q[15];
cx q[17], q[5];
cx q[17], q[3];
cx q[17], q[0];
cx q[18], q[15];
cx q[18], q[0];
cx q[19], q[15];
cx q[19], q[1];
cx q[19], q[0];
cx q[20], q[15];
cx q[20], q[5];
cx q[20], q[3];
cx q[21], q[15];
cx q[22], q[5];
cx q[22], q[1];
cx q[23], q[3];
cx q[24], q[5];
cx q[24], q[1];
cx q[24], q[0];
cx q[26], q[5];
cx q[28], q[5];
cx q[28], q[1];
cx q[29], q[1];
cx q[29], q[0];
cx q[27], q[3];
cx q[31], q[15];
cx q[31], q[0];
cx q[32], q[3];
cx q[32], q[0];
cx q[33], q[30];
cx q[33], q[15];
cx q[33], q[5];
cx q[33], q[3];
cx q[33], q[0];
cx q[34], q[15];
cx q[34], q[5];
cx q[34], q[3];
cx q[34], q[1];
cx q[35], q[5];
cx q[35], q[1];
cx q[22], q[30];
cx q[20], q[30];
cx q[19], q[30];
cx q[17], q[30];
cx q[16], q[30];
cx q[12], q[30];
cx q[12], q[15];
cx q[10], q[30];
cx q[10], q[15];
cx q[9], q[30];
cx q[8], q[30];
cx q[7], q[15];
cx q[6], q[30];
cx q[6], q[15];
cx q[4], q[30];
cx q[4], q[5];
cx q[3], q[15];
cx q[2], q[30];
cx q[2], q[5];
cx q[1], q[30];
cx q[1], q[15];
cx q[1], q[5];
