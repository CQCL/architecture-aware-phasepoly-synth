// Initial wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]
// Resulting wiring: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]
OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
cx q[9], q[32];
cx q[7], q[32];
cx q[25], q[32];
cx q[29], q[32];
cx q[37], q[32];
cx q[41], q[32];
cx q[0], q[32];
cx q[18], q[32];
cx q[30], q[32];
cx q[36], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[4], q[32];
cx q[5], q[32];
cx q[8], q[32];
cx q[19], q[32];
cx q[23], q[32];
cx q[27], q[32];
cx q[38], q[32];
cx q[42], q[32];
cx q[49], q[32];
cx q[51], q[32];
rz(0.25*pi) q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
cx q[4], q[32];
cx q[5], q[32];
cx q[8], q[32];
cx q[12], q[32];
cx q[13], q[32];
cx q[14], q[32];
cx q[19], q[32];
cx q[20], q[32];
cx q[21], q[32];
cx q[23], q[32];
cx q[26], q[32];
cx q[27], q[32];
cx q[28], q[32];
cx q[31], q[32];
cx q[33], q[32];
cx q[35], q[32];
cx q[38], q[32];
cx q[42], q[32];
cx q[45], q[32];
cx q[47], q[32];
cx q[49], q[32];
cx q[50], q[32];
cx q[51], q[32];
rz(0.25*pi) q[32];
cx q[0], q[32];
cx q[2], q[32];
cx q[5], q[32];
cx q[12], q[32];
cx q[13], q[32];
cx q[15], q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
cx q[21], q[32];
cx q[22], q[32];
cx q[23], q[32];
cx q[24], q[32];
cx q[26], q[32];
cx q[27], q[32];
cx q[30], q[32];
cx q[36], q[32];
cx q[43], q[32];
cx q[45], q[32];
cx q[51], q[32];
rz(0.25*pi) q[32];
cx q[7], q[32];
cx q[10], q[32];
cx q[2], q[32];
cx q[27], q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[3], q[32];
cx q[6], q[32];
cx q[11], q[32];
cx q[12], q[32];
cx q[13], q[32];
cx q[14], q[32];
cx q[15], q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[21], q[32];
cx q[22], q[32];
cx q[23], q[32];
cx q[24], q[32];
cx q[26], q[32];
cx q[29], q[32];
cx q[30], q[32];
cx q[31], q[32];
cx q[38], q[32];
cx q[39], q[32];
cx q[41], q[32];
cx q[42], q[32];
cx q[45], q[32];
cx q[47], q[32];
cx q[48], q[32];
cx q[49], q[32];
rz(0.25*pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[15], q[32];
cx q[23], q[32];
cx q[47], q[32];
cx q[1], q[32];
cx q[3], q[32];
cx q[8], q[32];
cx q[12], q[32];
cx q[13], q[32];
cx q[17], q[32];
cx q[20], q[32];
cx q[21], q[32];
cx q[22], q[32];
cx q[24], q[32];
cx q[25], q[32];
cx q[26], q[32];
cx q[28], q[32];
cx q[29], q[32];
cx q[36], q[32];
cx q[39], q[32];
cx q[40], q[32];
cx q[44], q[32];
cx q[46], q[32];
cx q[50], q[32];
cx q[52], q[32];
rz(0.25*pi) q[32];
cx q[1], q[32];
cx q[12], q[32];
cx q[14], q[32];
cx q[29], q[32];
cx q[51], q[32];
cx q[52], q[32];
cx q[3], q[32];
cx q[5], q[32];
cx q[11], q[32];
cx q[13], q[32];
cx q[16], q[32];
cx q[19], q[32];
cx q[21], q[32];
cx q[22], q[32];
cx q[24], q[32];
cx q[26], q[32];
cx q[28], q[32];
cx q[31], q[32];
cx q[33], q[32];
cx q[34], q[32];
cx q[37], q[32];
cx q[39], q[32];
cx q[40], q[32];
cx q[41], q[32];
cx q[43], q[32];
cx q[45], q[32];
cx q[48], q[32];
cx q[49], q[32];
rz(0.25*pi) q[32];
cx q[3], q[32];
cx q[8], q[32];
cx q[13], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
cx q[21], q[32];
cx q[26], q[32];
cx q[34], q[32];
cx q[35], q[32];
cx q[36], q[32];
cx q[42], q[32];
cx q[43], q[32];
cx q[50], q[32];
cx q[5], q[32];
cx q[11], q[32];
cx q[22], q[32];
cx q[31], q[32];
cx q[33], q[32];
cx q[37], q[32];
cx q[38], q[32];
cx q[39], q[32];
cx q[40], q[32];
cx q[41], q[32];
cx q[44], q[32];
cx q[46], q[32];
cx q[48], q[32];
rz(0.25*pi) q[32];
cx q[5], q[32];
cx q[6], q[32];
cx q[11], q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[22], q[32];
cx q[24], q[32];
cx q[25], q[32];
cx q[28], q[32];
cx q[30], q[32];
cx q[31], q[32];
cx q[33], q[32];
cx q[37], q[32];
cx q[38], q[32];
cx q[39], q[32];
cx q[40], q[32];
cx q[41], q[32];
cx q[44], q[32];
cx q[45], q[32];
cx q[46], q[32];
cx q[48], q[32];
cx q[49], q[32];
rz(0.25*pi) q[32];
cx q[1], q[32];
cx q[3], q[32];
cx q[4], q[32];
cx q[11], q[32];
cx q[12], q[32];
cx q[13], q[32];
cx q[14], q[32];
cx q[15], q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[22], q[32];
cx q[23], q[32];
cx q[24], q[32];
cx q[25], q[32];
cx q[29], q[32];
cx q[31], q[32];
cx q[34], q[32];
cx q[35], q[32];
cx q[43], q[32];
cx q[44], q[32];
cx q[45], q[32];
cx q[47], q[32];
cx q[48], q[32];
cx q[50], q[32];
cx q[51], q[32];
cx q[52], q[32];
rz(0.25*pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[5], q[32];
cx q[6], q[32];
cx q[8], q[32];
cx q[18], q[32];
cx q[19], q[32];
cx q[20], q[32];
cx q[22], q[32];
cx q[24], q[32];
cx q[26], q[32];
cx q[27], q[32];
cx q[28], q[32];
cx q[29], q[32];
cx q[31], q[32];
cx q[36], q[32];
cx q[41], q[32];
cx q[45], q[32];
cx q[48], q[32];
cx q[52], q[32];
rz(0.25*pi) q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[4], q[32];
cx q[11], q[32];
cx q[15], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[21], q[32];
cx q[22], q[32];
cx q[23], q[32];
cx q[24], q[32];
cx q[26], q[32];
cx q[28], q[32];
cx q[30], q[32];
cx q[34], q[32];
cx q[35], q[32];
cx q[37], q[32];
cx q[38], q[32];
cx q[42], q[32];
cx q[43], q[32];
cx q[46], q[32];
cx q[48], q[32];
cx q[50], q[32];
cx q[52], q[32];
rz(0.25*pi) q[32];
cx q[1], q[32];
cx q[4], q[32];
cx q[6], q[32];
cx q[8], q[32];
cx q[10], q[32];
cx q[11], q[32];
cx q[14], q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
cx q[23], q[32];
cx q[24], q[32];
cx q[26], q[32];
cx q[28], q[32];
cx q[29], q[32];
cx q[34], q[32];
cx q[35], q[32];
cx q[36], q[32];
cx q[37], q[32];
cx q[39], q[32];
cx q[40], q[32];
cx q[42], q[32];
cx q[45], q[32];
cx q[46], q[32];
cx q[47], q[32];
cx q[50], q[32];
cx q[51], q[32];
cx q[52], q[32];
rz(0.25*pi) q[32];
cx q[9], q[32];
cx q[48], q[32];
cx q[52], q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[3], q[32];
cx q[4], q[32];
cx q[5], q[32];
cx q[7], q[32];
cx q[10], q[32];
cx q[11], q[32];
cx q[13], q[32];
cx q[15], q[32];
cx q[16], q[32];
cx q[18], q[32];
cx q[19], q[32];
cx q[21], q[32];
cx q[22], q[32];
cx q[24], q[32];
cx q[26], q[32];
cx q[27], q[32];
cx q[29], q[32];
cx q[30], q[32];
cx q[38], q[32];
cx q[40], q[32];
cx q[42], q[32];
cx q[45], q[32];
cx q[46], q[32];
cx q[47], q[32];
cx q[51], q[32];
rz(0.25*pi) q[32];
cx q[0], q[32];
cx q[6], q[32];
cx q[24], q[32];
cx q[3], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[10], q[32];
cx q[11], q[32];
cx q[15], q[32];
cx q[17], q[32];
cx q[21], q[32];
cx q[22], q[32];
cx q[25], q[32];
cx q[27], q[32];
cx q[29], q[32];
cx q[30], q[32];
cx q[34], q[32];
cx q[35], q[32];
cx q[36], q[32];
cx q[41], q[32];
cx q[42], q[32];
cx q[43], q[32];
cx q[45], q[32];
cx q[46], q[32];
cx q[49], q[32];
rz(0.25*pi) q[32];
cx q[4], q[32];
cx q[13], q[32];
cx q[21], q[32];
cx q[22], q[32];
cx q[23], q[32];
cx q[26], q[32];
cx q[29], q[32];
cx q[33], q[32];
cx q[41], q[32];
cx q[45], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[7], q[32];
cx q[8], q[32];
cx q[16], q[32];
cx q[19], q[32];
cx q[42], q[32];
cx q[49], q[32];
cx q[50], q[32];
cx q[3], q[32];
cx q[5], q[32];
cx q[10], q[32];
cx q[11], q[32];
cx q[12], q[32];
cx q[15], q[32];
cx q[17], q[32];
cx q[37], q[32];
cx q[39], q[32];
rz(0.25*pi) q[32];
cx q[3], q[32];
cx q[5], q[32];
cx q[10], q[32];
cx q[11], q[32];
cx q[12], q[32];
cx q[14], q[32];
cx q[15], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[20], q[32];
cx q[25], q[32];
cx q[27], q[32];
cx q[30], q[32];
cx q[31], q[32];
cx q[35], q[32];
cx q[37], q[32];
cx q[39], q[32];
cx q[40], q[32];
cx q[43], q[32];
cx q[44], q[32];
cx q[46], q[32];
rz(0.25*pi) q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
cx q[5], q[32];
cx q[7], q[32];
cx q[8], q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[19], q[32];
cx q[25], q[32];
cx q[30], q[32];
cx q[34], q[32];
cx q[35], q[32];
cx q[36], q[32];
cx q[37], q[32];
cx q[38], q[32];
cx q[39], q[32];
cx q[42], q[32];
cx q[47], q[32];
cx q[49], q[32];
cx q[50], q[32];
cx q[51], q[32];
rz(0.25*pi) q[32];
cx q[0], q[6];
cx q[2], q[6];
cx q[3], q[6];
cx q[4], q[6];
cx q[14], q[6];
cx q[15], q[6];
cx q[20], q[6];
cx q[22], q[6];
cx q[23], q[6];
cx q[31], q[6];
cx q[33], q[6];
cx q[40], q[6];
cx q[42], q[6];
cx q[46], q[6];
cx q[48], q[6];
cx q[1], q[6];
cx q[11], q[6];
cx q[21], q[6];
cx q[34], q[6];
cx q[35], q[6];
cx q[38], q[6];
cx q[50], q[6];
rz(0.25*pi) q[6];
cx q[1], q[6];
cx q[5], q[6];
cx q[7], q[6];
cx q[9], q[6];
cx q[10], q[6];
cx q[11], q[6];
cx q[16], q[6];
cx q[17], q[6];
cx q[21], q[6];
cx q[25], q[6];
cx q[26], q[6];
cx q[28], q[6];
cx q[34], q[6];
cx q[35], q[6];
cx q[36], q[6];
cx q[37], q[6];
cx q[38], q[6];
cx q[39], q[6];
cx q[43], q[6];
cx q[44], q[6];
cx q[45], q[6];
cx q[49], q[6];
cx q[50], q[6];
cx q[52], q[6];
rz(0.25*pi) q[6];
cx q[4], q[6];
cx q[2], q[6];
cx q[21], q[6];
cx q[9], q[6];
cx q[10], q[6];
cx q[12], q[6];
cx q[16], q[6];
cx q[43], q[6];
cx q[48], q[6];
cx q[1], q[6];
cx q[5], q[6];
cx q[11], q[6];
cx q[14], q[6];
cx q[15], q[6];
cx q[18], q[6];
cx q[19], q[6];
cx q[20], q[6];
cx q[23], q[6];
cx q[28], q[6];
cx q[30], q[6];
cx q[31], q[6];
cx q[37], q[6];
cx q[38], q[6];
cx q[39], q[6];
cx q[44], q[6];
cx q[50], q[6];
rz(0.25*pi) q[6];
cx q[1], q[6];
cx q[5], q[6];
cx q[8], q[6];
cx q[14], q[6];
cx q[17], q[6];
cx q[18], q[6];
cx q[22], q[6];
cx q[24], q[6];
cx q[27], q[6];
cx q[33], q[6];
cx q[36], q[6];
cx q[37], q[6];
cx q[41], q[6];
cx q[44], q[6];
cx q[46], q[6];
cx q[51], q[6];
cx q[3], q[6];
cx q[7], q[6];
cx q[15], q[6];
cx q[19], q[6];
cx q[23], q[6];
cx q[28], q[6];
cx q[34], q[6];
cx q[38], q[6];
cx q[45], q[6];
rz(0.25*pi) q[6];
cx q[3], q[6];
cx q[7], q[6];
cx q[11], q[6];
cx q[15], q[6];
cx q[19], q[6];
cx q[20], q[6];
cx q[23], q[6];
cx q[25], q[6];
cx q[26], q[6];
cx q[28], q[6];
cx q[29], q[6];
cx q[30], q[6];
cx q[31], q[6];
cx q[34], q[6];
cx q[38], q[6];
cx q[39], q[6];
cx q[40], q[6];
cx q[45], q[6];
cx q[49], q[6];
cx q[50], q[6];
rz(0.25*pi) q[6];
cx q[1], q[6];
cx q[3], q[6];
cx q[5], q[6];
cx q[7], q[6];
cx q[8], q[6];
cx q[10], q[6];
cx q[12], q[6];
cx q[13], q[6];
cx q[14], q[6];
cx q[16], q[6];
cx q[17], q[6];
cx q[19], q[6];
cx q[20], q[6];
cx q[22], q[6];
cx q[24], q[6];
cx q[28], q[6];
cx q[30], q[6];
cx q[33], q[6];
cx q[35], q[6];
cx q[41], q[6];
cx q[42], q[6];
cx q[43], q[6];
cx q[44], q[6];
cx q[46], q[6];
cx q[47], q[6];
cx q[48], q[6];
cx q[49], q[6];
cx q[51], q[6];
rz(0.25*pi) q[6];
cx q[3], q[6];
cx q[5], q[6];
cx q[7], q[6];
cx q[8], q[6];
cx q[9], q[6];
cx q[10], q[6];
cx q[11], q[6];
cx q[12], q[6];
cx q[14], q[6];
cx q[16], q[6];
cx q[18], q[6];
cx q[23], q[6];
cx q[24], q[6];
cx q[26], q[6];
cx q[27], q[6];
cx q[28], q[6];
cx q[29], q[6];
cx q[30], q[6];
cx q[36], q[6];
cx q[37], q[6];
cx q[38], q[6];
cx q[39], q[6];
cx q[40], q[6];
cx q[43], q[6];
cx q[48], q[6];
cx q[52], q[6];
rz(0.25*pi) q[6];
cx q[1], q[6];
cx q[2], q[6];
cx q[3], q[6];
cx q[10], q[6];
cx q[11], q[6];
cx q[12], q[6];
cx q[14], q[6];
cx q[15], q[6];
cx q[18], q[6];
cx q[19], q[6];
cx q[21], q[6];
cx q[22], q[6];
cx q[23], q[6];
cx q[25], q[6];
cx q[27], q[6];
cx q[29], q[6];
cx q[31], q[6];
cx q[34], q[6];
cx q[36], q[6];
cx q[37], q[6];
cx q[40], q[6];
cx q[41], q[6];
cx q[42], q[6];
cx q[43], q[6];
cx q[45], q[6];
cx q[48], q[6];
cx q[49], q[6];
rz(0.25*pi) q[6];
cx q[0], q[6];
cx q[3], q[6];
cx q[5], q[6];
cx q[7], q[6];
cx q[19], q[6];
cx q[21], q[6];
cx q[25], q[6];
cx q[29], q[6];
cx q[31], q[6];
cx q[38], q[6];
cx q[42], q[6];
cx q[49], q[6];
cx q[1], q[6];
cx q[12], q[6];
cx q[14], q[6];
cx q[16], q[6];
cx q[17], q[6];
cx q[20], q[6];
cx q[22], q[6];
cx q[23], q[6];
cx q[28], q[6];
cx q[34], q[6];
cx q[35], q[6];
cx q[36], q[6];
cx q[37], q[6];
cx q[40], q[6];
cx q[41], q[6];
cx q[44], q[6];
cx q[47], q[6];
rz(0.25*pi) q[6];
cx q[1], q[6];
cx q[10], q[6];
cx q[11], q[6];
cx q[12], q[6];
cx q[14], q[6];
cx q[15], q[6];
cx q[16], q[6];
cx q[17], q[6];
cx q[20], q[6];
cx q[22], q[6];
cx q[23], q[6];
cx q[24], q[6];
cx q[27], q[6];
cx q[28], q[6];
cx q[30], q[6];
cx q[33], q[6];
cx q[34], q[6];
cx q[35], q[6];
cx q[36], q[6];
cx q[37], q[6];
cx q[40], q[6];
cx q[41], q[6];
cx q[43], q[6];
cx q[44], q[6];
cx q[46], q[6];
cx q[47], q[6];
cx q[48], q[6];
cx q[50], q[6];
rz(0.25*pi) q[6];
cx q[52], q[11];
cx q[0], q[11];
cx q[38], q[11];
cx q[1], q[11];
cx q[2], q[11];
cx q[4], q[11];
cx q[5], q[11];
cx q[10], q[11];
cx q[13], q[11];
cx q[16], q[11];
cx q[18], q[11];
cx q[21], q[11];
cx q[25], q[11];
cx q[26], q[11];
cx q[28], q[11];
cx q[30], q[11];
cx q[31], q[11];
cx q[33], q[11];
cx q[34], q[11];
cx q[36], q[11];
cx q[39], q[11];
cx q[41], q[11];
cx q[43], q[11];
cx q[44], q[11];
cx q[45], q[11];
cx q[46], q[11];
cx q[47], q[11];
cx q[50], q[11];
cx q[51], q[11];
rz(0.25*pi) q[11];
cx q[1], q[11];
cx q[3], q[11];
cx q[10], q[11];
cx q[17], q[11];
cx q[30], q[11];
cx q[2], q[11];
cx q[19], q[11];
cx q[22], q[11];
cx q[4], q[11];
cx q[9], q[11];
cx q[20], q[11];
cx q[33], q[11];
cx q[34], q[11];
cx q[37], q[11];
cx q[48], q[11];
cx q[5], q[11];
cx q[12], q[11];
cx q[14], q[11];
cx q[21], q[11];
cx q[26], q[11];
cx q[31], q[11];
cx q[39], q[11];
cx q[40], q[11];
cx q[41], q[11];
cx q[46], q[11];
cx q[47], q[11];
cx q[50], q[11];
rz(0.25*pi) q[11];
cx q[5], q[11];
cx q[8], q[11];
cx q[12], q[11];
cx q[13], q[11];
cx q[14], q[11];
cx q[16], q[11];
cx q[18], q[11];
cx q[21], q[11];
cx q[25], q[11];
cx q[26], q[11];
cx q[28], q[11];
cx q[29], q[11];
cx q[31], q[11];
cx q[35], q[11];
cx q[36], q[11];
cx q[39], q[11];
cx q[40], q[11];
cx q[41], q[11];
cx q[43], q[11];
cx q[46], q[11];
cx q[47], q[11];
cx q[50], q[11];
cx q[51], q[11];
rz(0.25*pi) q[11];
cx q[4], q[11];
cx q[7], q[11];
cx q[9], q[11];
cx q[14], q[11];
cx q[15], q[11];
cx q[16], q[11];
cx q[20], q[11];
cx q[24], q[11];
cx q[25], q[11];
cx q[27], q[11];
cx q[28], q[11];
cx q[29], q[11];
cx q[33], q[11];
cx q[34], q[11];
cx q[36], q[11];
cx q[37], q[11];
cx q[39], q[11];
cx q[40], q[11];
cx q[41], q[11];
cx q[45], q[11];
cx q[46], q[11];
cx q[48], q[11];
cx q[49], q[11];
cx q[50], q[11];
rz(0.25*pi) q[11];
cx q[2], q[11];
cx q[4], q[11];
cx q[7], q[11];
cx q[9], q[11];
cx q[12], q[11];
cx q[13], q[11];
cx q[14], q[11];
cx q[15], q[11];
cx q[18], q[11];
cx q[19], q[11];
cx q[20], q[11];
cx q[22], q[11];
cx q[23], q[11];
cx q[24], q[11];
cx q[25], q[11];
cx q[26], q[11];
cx q[27], q[11];
cx q[29], q[11];
cx q[33], q[11];
cx q[37], q[11];
cx q[39], q[11];
cx q[42], q[11];
cx q[44], q[11];
cx q[47], q[11];
cx q[48], q[11];
rz(0.25*pi) q[11];
cx q[0], q[11];
cx q[1], q[11];
cx q[9], q[11];
cx q[10], q[11];
cx q[12], q[11];
cx q[14], q[11];
cx q[15], q[11];
cx q[18], q[11];
cx q[19], q[11];
cx q[20], q[11];
cx q[22], q[11];
cx q[23], q[11];
cx q[25], q[11];
cx q[29], q[11];
cx q[31], q[11];
cx q[34], q[11];
cx q[36], q[11];
cx q[38], q[11];
cx q[41], q[11];
cx q[43], q[11];
cx q[45], q[11];
cx q[46], q[11];
cx q[47], q[11];
cx q[49], q[11];
rz(0.25*pi) q[11];
cx q[8], q[14];
cx q[16], q[14];
cx q[20], q[14];
cx q[52], q[14];
cx q[0], q[14];
cx q[1], q[14];
cx q[7], q[14];
cx q[10], q[14];
cx q[12], q[14];
cx q[13], q[14];
cx q[17], q[14];
cx q[18], q[14];
cx q[19], q[14];
cx q[21], q[14];
cx q[27], q[14];
cx q[29], q[14];
cx q[30], q[14];
cx q[31], q[14];
cx q[38], q[14];
cx q[44], q[14];
cx q[46], q[14];
cx q[47], q[14];
rz(0.25*pi) q[14];
cx q[1], q[14];
cx q[2], q[14];
cx q[10], q[14];
cx q[0], q[14];
cx q[21], q[14];
cx q[27], q[14];
cx q[45], q[14];
cx q[3], q[14];
cx q[4], q[14];
cx q[5], q[14];
cx q[7], q[14];
cx q[12], q[14];
cx q[13], q[14];
cx q[15], q[14];
cx q[17], q[14];
cx q[22], q[14];
cx q[23], q[14];
cx q[24], q[14];
cx q[26], q[14];
cx q[35], q[14];
cx q[38], q[14];
cx q[41], q[14];
cx q[42], q[14];
cx q[46], q[14];
cx q[49], q[14];
cx q[50], q[14];
cx q[51], q[14];
rz(0.25*pi) q[14];
cx q[3], q[14];
cx q[4], q[14];
cx q[5], q[14];
cx q[7], q[14];
cx q[9], q[14];
cx q[12], q[14];
cx q[13], q[14];
cx q[15], q[14];
cx q[17], q[14];
cx q[18], q[14];
cx q[22], q[14];
cx q[23], q[14];
cx q[24], q[14];
cx q[26], q[14];
cx q[28], q[14];
cx q[29], q[14];
cx q[31], q[14];
cx q[34], q[14];
cx q[35], q[14];
cx q[36], q[14];
cx q[37], q[14];
cx q[38], q[14];
cx q[41], q[14];
cx q[42], q[14];
cx q[43], q[14];
cx q[44], q[14];
cx q[46], q[14];
cx q[47], q[14];
cx q[48], q[14];
cx q[49], q[14];
cx q[50], q[14];
cx q[51], q[14];
rz(0.25*pi) q[14];
cx q[0], q[14];
cx q[3], q[14];
cx q[12], q[14];
cx q[13], q[14];
cx q[17], q[14];
cx q[18], q[14];
cx q[19], q[14];
cx q[21], q[14];
cx q[22], q[14];
cx q[24], q[14];
cx q[27], q[14];
cx q[29], q[14];
cx q[30], q[14];
cx q[34], q[14];
cx q[35], q[14];
cx q[37], q[14];
cx q[39], q[14];
cx q[40], q[14];
cx q[41], q[14];
cx q[42], q[14];
cx q[43], q[14];
cx q[44], q[14];
cx q[45], q[14];
cx q[46], q[14];
cx q[47], q[14];
cx q[51], q[14];
rz(0.25*pi) q[14];
cx q[13], q[3];
cx q[0], q[3];
cx q[16], q[3];
cx q[31], q[3];
cx q[37], q[3];
cx q[51], q[3];
cx q[1], q[3];
cx q[10], q[3];
cx q[15], q[3];
cx q[28], q[3];
cx q[29], q[3];
cx q[38], q[3];
cx q[48], q[3];
cx q[2], q[3];
cx q[5], q[3];
cx q[17], q[3];
cx q[19], q[3];
cx q[20], q[3];
cx q[25], q[3];
cx q[26], q[3];
cx q[27], q[3];
cx q[34], q[3];
cx q[35], q[3];
cx q[39], q[3];
cx q[41], q[3];
cx q[44], q[3];
cx q[45], q[3];
cx q[46], q[3];
cx q[47], q[3];
cx q[49], q[3];
cx q[50], q[3];
rz(0.25*pi) q[3];
cx q[2], q[3];
cx q[4], q[3];
cx q[5], q[3];
cx q[9], q[3];
cx q[12], q[3];
cx q[17], q[3];
cx q[19], q[3];
cx q[20], q[3];
cx q[21], q[3];
cx q[23], q[3];
cx q[25], q[3];
cx q[26], q[3];
cx q[27], q[3];
cx q[30], q[3];
cx q[34], q[3];
cx q[35], q[3];
cx q[36], q[3];
cx q[39], q[3];
cx q[41], q[3];
cx q[43], q[3];
cx q[44], q[3];
cx q[45], q[3];
cx q[46], q[3];
cx q[47], q[3];
cx q[49], q[3];
cx q[50], q[3];
rz(0.25*pi) q[3];
cx q[1], q[3];
cx q[9], q[3];
cx q[10], q[3];
cx q[12], q[3];
cx q[15], q[3];
cx q[18], q[3];
cx q[19], q[3];
cx q[20], q[3];
cx q[21], q[3];
cx q[22], q[3];
cx q[28], q[3];
cx q[29], q[3];
cx q[38], q[3];
cx q[40], q[3];
cx q[41], q[3];
cx q[42], q[3];
cx q[44], q[3];
cx q[45], q[3];
cx q[48], q[3];
cx q[49], q[3];
cx q[50], q[3];
cx q[52], q[3];
rz(0.25*pi) q[3];
cx q[0], q[3];
cx q[1], q[3];
cx q[2], q[3];
cx q[5], q[3];
cx q[8], q[3];
cx q[15], q[3];
cx q[16], q[3];
cx q[18], q[3];
cx q[19], q[3];
cx q[21], q[3];
cx q[23], q[3];
cx q[24], q[3];
cx q[25], q[3];
cx q[27], q[3];
cx q[30], q[3];
cx q[31], q[3];
cx q[33], q[3];
cx q[35], q[3];
cx q[36], q[3];
cx q[37], q[3];
cx q[38], q[3];
cx q[40], q[3];
cx q[41], q[3];
cx q[47], q[3];
cx q[50], q[3];
cx q[51], q[3];
rz(0.25*pi) q[3];
cx q[7], q[5];
cx q[13], q[5];
cx q[18], q[5];
cx q[21], q[5];
cx q[29], q[5];
cx q[30], q[5];
cx q[39], q[5];
cx q[42], q[5];
cx q[45], q[5];
cx q[0], q[5];
cx q[4], q[5];
cx q[10], q[5];
cx q[16], q[5];
cx q[19], q[5];
cx q[22], q[5];
cx q[26], q[5];
cx q[28], q[5];
cx q[31], q[5];
cx q[33], q[5];
cx q[35], q[5];
cx q[44], q[5];
cx q[48], q[5];
rz(0.25*pi) q[5];
cx q[0], q[5];
cx q[1], q[5];
cx q[2], q[5];
cx q[4], q[5];
cx q[9], q[5];
cx q[10], q[5];
cx q[15], q[5];
cx q[16], q[5];
cx q[17], q[5];
cx q[19], q[5];
cx q[20], q[5];
cx q[22], q[5];
cx q[23], q[5];
cx q[24], q[5];
cx q[25], q[5];
cx q[26], q[5];
cx q[27], q[5];
cx q[28], q[5];
cx q[31], q[5];
cx q[33], q[5];
cx q[35], q[5];
cx q[37], q[5];
cx q[40], q[5];
cx q[44], q[5];
cx q[46], q[5];
cx q[47], q[5];
cx q[48], q[5];
cx q[49], q[5];
cx q[50], q[5];
rz(0.25*pi) q[5];
cx q[1], q[2];
cx q[7], q[2];
cx q[24], q[2];
cx q[25], q[2];
cx q[27], q[2];
cx q[40], q[2];
cx q[43], q[2];
cx q[44], q[2];
cx q[46], q[2];
cx q[52], q[2];
cx q[8], q[2];
cx q[10], q[2];
cx q[15], q[2];
cx q[22], q[2];
cx q[23], q[2];
cx q[38], q[2];
cx q[41], q[2];
cx q[49], q[2];
rz(0.25*pi) q[2];
cx q[8], q[2];
cx q[10], q[2];
cx q[12], q[2];
cx q[13], q[2];
cx q[15], q[2];
cx q[16], q[2];
cx q[17], q[2];
cx q[19], q[2];
cx q[20], q[2];
cx q[21], q[2];
cx q[22], q[2];
cx q[23], q[2];
cx q[26], q[2];
cx q[28], q[2];
cx q[30], q[2];
cx q[31], q[2];
cx q[34], q[2];
cx q[36], q[2];
cx q[37], q[2];
cx q[38], q[2];
cx q[39], q[2];
cx q[41], q[2];
cx q[45], q[2];
cx q[48], q[2];
cx q[49], q[2];
cx q[50], q[2];
cx q[51], q[2];
rz(0.25*pi) q[2];
cx q[0], q[1];
cx q[4], q[1];
cx q[8], q[1];
cx q[9], q[1];
cx q[12], q[1];
cx q[15], q[1];
cx q[18], q[1];
cx q[19], q[1];
cx q[20], q[1];
cx q[23], q[1];
cx q[28], q[1];
cx q[29], q[1];
cx q[33], q[1];
cx q[35], q[1];
cx q[37], q[1];
cx q[39], q[1];
cx q[41], q[1];
cx q[42], q[1];
cx q[46], q[1];
cx q[47], q[1];
cx q[48], q[1];
cx q[51], q[1];
rz(0.25*pi) q[1];
cx q[15], q[4];
cx q[28], q[4];
cx q[33], q[4];
cx q[38], q[4];
cx q[0], q[4];
cx q[8], q[4];
cx q[13], q[4];
cx q[16], q[4];
cx q[18], q[4];
cx q[19], q[4];
cx q[21], q[4];
cx q[23], q[4];
cx q[26], q[4];
cx q[35], q[4];
cx q[39], q[4];
cx q[42], q[4];
cx q[46], q[4];
cx q[47], q[4];
cx q[50], q[4];
cx q[51], q[4];
rz(0.25*pi) q[4];
cx q[0], q[4];
cx q[20], q[4];
cx q[21], q[4];
cx q[23], q[4];
cx q[27], q[4];
cx q[30], q[4];
cx q[44], q[4];
cx q[45], q[4];
cx q[46], q[4];
cx q[8], q[4];
cx q[9], q[4];
cx q[12], q[4];
cx q[16], q[4];
cx q[18], q[4];
cx q[19], q[4];
cx q[22], q[4];
cx q[25], q[4];
cx q[26], q[4];
cx q[34], q[4];
cx q[41], q[4];
cx q[42], q[4];
cx q[43], q[4];
cx q[48], q[4];
cx q[49], q[4];
cx q[50], q[4];
rz(0.25*pi) q[4];
cx q[8], q[4];
cx q[9], q[4];
cx q[12], q[4];
cx q[13], q[4];
cx q[16], q[4];
cx q[17], q[4];
cx q[18], q[4];
cx q[19], q[4];
cx q[22], q[4];
cx q[25], q[4];
cx q[26], q[4];
cx q[29], q[4];
cx q[34], q[4];
cx q[35], q[4];
cx q[36], q[4];
cx q[39], q[4];
cx q[41], q[4];
cx q[42], q[4];
cx q[43], q[4];
cx q[47], q[4];
cx q[48], q[4];
cx q[49], q[4];
cx q[50], q[4];
cx q[51], q[4];
rz(0.25*pi) q[4];
cx q[7], q[0];
cx q[10], q[0];
cx q[12], q[0];
cx q[17], q[0];
cx q[18], q[0];
cx q[19], q[0];
cx q[20], q[0];
cx q[21], q[0];
cx q[23], q[0];
cx q[25], q[0];
cx q[29], q[0];
cx q[30], q[0];
cx q[31], q[0];
cx q[34], q[0];
cx q[38], q[0];
cx q[39], q[0];
cx q[41], q[0];
cx q[43], q[0];
cx q[45], q[0];
cx q[48], q[0];
cx q[50], q[0];
rz(0.25*pi) q[0];
cx q[2], q[0];
cx q[4], q[3];
cx q[4], q[1];
cx q[5], q[3];
cx q[7], q[6];
cx q[7], q[5];
cx q[7], q[3];
cx q[7], q[2];
cx q[8], q[6];
cx q[8], q[4];
cx q[9], q[5];
cx q[9], q[3];
cx q[10], q[6];
cx q[10], q[0];
cx q[11], q[6];
cx q[11], q[5];
cx q[11], q[4];
cx q[11], q[2];
cx q[11], q[0];
cx q[9], q[1];
cx q[12], q[2];
cx q[12], q[1];
cx q[13], q[11];
cx q[13], q[4];
cx q[14], q[11];
cx q[14], q[6];
cx q[14], q[5];
cx q[14], q[4];
cx q[14], q[2];
cx q[14], q[0];
cx q[15], q[11];
cx q[15], q[6];
cx q[15], q[3];
cx q[15], q[2];
cx q[15], q[0];
cx q[16], q[14];
cx q[16], q[4];
cx q[16], q[3];
cx q[16], q[2];
cx q[16], q[1];
cx q[16], q[0];
cx q[17], q[11];
cx q[17], q[1];
cx q[17], q[0];
cx q[18], q[14];
cx q[18], q[11];
cx q[18], q[6];
cx q[18], q[2];
cx q[19], q[11];
cx q[19], q[6];
cx q[19], q[5];
cx q[19], q[3];
cx q[19], q[0];
cx q[20], q[14];
cx q[20], q[11];
cx q[20], q[6];
cx q[20], q[3];
cx q[20], q[0];
cx q[21], q[14];
cx q[21], q[6];
cx q[21], q[5];
cx q[21], q[2];
cx q[22], q[14];
cx q[22], q[3];
cx q[23], q[5];
cx q[23], q[3];
cx q[23], q[0];
cx q[23], q[1];
cx q[11], q[4];
cx q[24], q[14];
cx q[24], q[11];
cx q[24], q[6];
cx q[25], q[11];
cx q[25], q[6];
cx q[25], q[0];
cx q[26], q[6];
cx q[26], q[4];
cx q[26], q[3];
cx q[26], q[2];
cx q[26], q[1];
cx q[26], q[0];
cx q[27], q[14];
cx q[27], q[11];
cx q[27], q[6];
cx q[27], q[4];
cx q[27], q[3];
cx q[27], q[1];
cx q[28], q[14];
cx q[28], q[4];
cx q[28], q[3];
cx q[28], q[2];
cx q[28], q[0];
cx q[29], q[14];
cx q[29], q[11];
cx q[29], q[6];
cx q[29], q[4];
cx q[29], q[2];
cx q[11], q[5];
cx q[30], q[6];
cx q[30], q[5];
cx q[30], q[4];
cx q[30], q[2];
cx q[30], q[1];
cx q[31], q[6];
cx q[31], q[2];
cx q[32], q[5];
cx q[33], q[32];
cx q[33], q[5];
cx q[33], q[4];
cx q[34], q[32];
cx q[34], q[6];
cx q[34], q[5];
cx q[34], q[2];
cx q[35], q[6];
cx q[35], q[1];
cx q[2], q[0];
cx q[35], q[3];
cx q[35], q[14];
cx q[36], q[14];
cx q[36], q[11];
cx q[36], q[6];
cx q[36], q[4];
cx q[36], q[3];
cx q[36], q[1];
cx q[37], q[32];
cx q[37], q[11];
cx q[37], q[3];
cx q[37], q[1];
cx q[38], q[14];
cx q[38], q[11];
cx q[38], q[4];
cx q[38], q[2];
cx q[38], q[1];
cx q[38], q[0];
cx q[39], q[32];
cx q[39], q[14];
cx q[39], q[3];
cx q[39], q[2];
cx q[39], q[1];
cx q[39], q[0];
cx q[40], q[32];
cx q[40], q[14];
cx q[40], q[3];
cx q[40], q[2];
cx q[41], q[14];
cx q[41], q[6];
cx q[41], q[1];
cx q[41], q[0];
cx q[32], q[5];
cx q[42], q[14];
cx q[42], q[5];
cx q[42], q[4];
cx q[42], q[3];
cx q[43], q[32];
cx q[43], q[11];
cx q[43], q[6];
cx q[43], q[3];
cx q[43], q[0];
cx q[44], q[32];
cx q[44], q[14];
cx q[44], q[11];
cx q[44], q[6];
cx q[44], q[4];
cx q[44], q[1];
cx q[45], q[32];
cx q[45], q[11];
cx q[45], q[6];
cx q[45], q[5];
cx q[45], q[4];
cx q[45], q[3];
cx q[45], q[1];
cx q[45], q[0];
cx q[46], q[11];
cx q[46], q[6];
cx q[46], q[5];
cx q[46], q[3];
cx q[46], q[1];
cx q[47], q[32];
cx q[47], q[14];
cx q[47], q[6];
cx q[47], q[5];
cx q[47], q[1];
cx q[11], q[2];
cx q[48], q[14];
cx q[48], q[2];
cx q[48], q[1];
cx q[48], q[0];
cx q[49], q[6];
cx q[49], q[5];
cx q[50], q[6];
cx q[50], q[5];
cx q[50], q[2];
cx q[50], q[1];
cx q[50], q[0];
cx q[51], q[14];
cx q[51], q[6];
cx q[51], q[2];
cx q[51], q[1];
cx q[52], q[14];
cx q[52], q[2];
cx q[52], q[3];
cx q[50], q[4];
cx q[51], q[11];
cx q[52], q[32];
cx q[30], q[32];
cx q[28], q[32];
cx q[26], q[32];
cx q[25], q[32];
cx q[24], q[32];
cx q[23], q[32];
cx q[20], q[32];
cx q[17], q[32];
cx q[16], q[32];
cx q[15], q[32];
cx q[13], q[32];
cx q[11], q[32];
cx q[10], q[32];
cx q[10], q[11];
cx q[9], q[32];
cx q[9], q[14];
cx q[9], q[11];
cx q[8], q[32];
cx q[8], q[14];
cx q[7], q[14];
cx q[7], q[11];
cx q[6], q[32];
cx q[5], q[32];
cx q[5], q[11];
cx q[5], q[6];
cx q[4], q[32];
cx q[3], q[14];
cx q[3], q[6];
cx q[2], q[32];
cx q[2], q[14];
cx q[2], q[11];
cx q[2], q[5];
cx q[1], q[14];
cx q[1], q[11];
cx q[1], q[2];
cx q[0], q[32];
cx q[0], q[6];
cx q[0], q[2];
cx q[0], q[1];
