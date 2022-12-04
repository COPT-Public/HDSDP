int nCol = 227; 
int nRow = 146; 
 int Ap[] = {0, 5, 10, 11, 18, 23, 28, 34, 37, 40, 43, 47, 51, 54, 57, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 126, 129, 133, 136, 140, 144, 148, 152, 156, 160, 164, 168, 175, 182, 186, 190, 194, 202, 209, 216, 220, 224, 228, 231, 236, 241, 250, 259, 268, 280, 292, 296, 300, 303, 306, 313, 320, 325, 330, 335, 342, 345, 354, 363, 372, 384, 396, 401, 406, 410, 414, 421, 425, 428, 431, 434, 437, 440, 443, 449, 456, 459, 462, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598}; 
int Ai[] = {13, 15, 16, 51, 53, 13, 15, 16, 51, 54, 13, 10, 16, 37, 49, 50, 51, 55, 6, 38, 39, 48, 56, 6, 39, 45, 48, 57, 10, 17, 37, 49, 50, 58, 17, 41, 59, 11, 17, 60, 2, 7, 61, 3, 8, 52, 62, 4, 9, 52, 63, 2, 18, 64, 3, 18, 65, 4, 18, 66, 2, 7, 19, 67, 3, 8, 19, 68, 4, 9, 19, 69, 2, 7, 20, 70, 3, 8, 20, 71, 4, 9, 20, 72, 2, 7, 22, 73, 3, 8, 22, 74, 4, 9, 22, 75, 2, 7, 23, 76, 3, 8, 23, 77, 4, 9, 23, 78, 2, 7, 24, 79, 3, 8, 24, 80, 4, 9, 24, 81, 10, 25, 37, 49, 50, 82, 25, 41, 83, 2, 7, 28, 84, 11, 25, 85, 3, 8, 28, 86, 4, 9, 28, 87, 2, 7, 26, 88, 3, 8, 26, 89, 4, 9, 26, 90, 2, 7, 27, 91, 3, 8, 27, 92, 4, 9, 27, 93, 0, 1, 12, 14, 21, 38, 94, 0, 1, 12, 14, 21, 45, 95, 2, 7, 29, 96, 3, 8, 29, 97, 4, 9, 29, 98, 1, 11, 12, 13, 17, 25, 31, 99, 1, 11, 12, 13, 25, 31, 100, 1, 11, 12, 13, 25, 31, 101, 1, 11, 12, 102, 1, 11, 12, 103, 11, 14, 32, 104, 32, 40, 105, 13, 25, 36, 51, 106, 25, 36, 46, 51, 107, 0, 1, 11, 12, 14, 21, 33, 38, 108, 0, 1, 11, 12, 14, 21, 33, 45, 109, 1, 2, 11, 12, 30, 38, 39, 51, 110, 1, 3, 11, 12, 13, 17, 30, 35, 38, 39, 51, 111, 1, 4, 11, 12, 13, 17, 34, 35, 38, 39, 51, 112, 1, 12, 38, 113, 1, 12, 39, 114, 39, 47, 115, 11, 39, 116, 1, 5, 12, 38, 39, 40, 117, 1, 5, 12, 39, 40, 45, 118, 17, 39, 42, 51, 119, 39, 42, 45, 51, 120, 38, 39, 42, 51, 121, 1, 11, 12, 33, 38, 43, 122, 43, 45, 123, 1, 6, 11, 12, 33, 38, 39, 48, 124, 1, 6, 11, 12, 33, 39, 45, 48, 125, 1, 7, 11, 12, 30, 39, 45, 51, 126, 1, 8, 11, 12, 13, 17, 30, 35, 39, 45, 51, 127, 1, 9, 11, 12, 13, 17, 34, 35, 39, 45, 51, 128, 16, 44, 46, 51, 129, 16, 44, 46, 51, 130, 16, 44, 51, 131, 1, 12, 45, 132, 10, 37, 46, 49, 50, 51, 133, 0, 14, 21, 134, 1, 12, 135, 1, 12, 136, 5, 40, 137, 6, 48, 138, 10, 49, 139, 10, 50, 140, 10, 37, 49, 50, 51, 141, 1, 11, 12, 16, 46, 51, 142, 41, 51, 143, 47, 51, 144, 11, 51, 145, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145}; 
double Ax[] = {0.506, 1.0, 0.182, 0.312, 1.0, 0.638, 1.0, 0.05, 0.312, 1.0, -1.0, -9.5, 0.92, 1.0, -0.042, -0.063, 0.08, 1.0, 1.0, 0.825, 0.175, -16.0, 1.0, 1.0, 0.175, 0.825, -21.0, 1.0, 3.6, 1.0, 1.0, -0.042, -0.063, 1.0, 1.0, 1.0, 1.0, -0.8, 1.0, 1.0, -1.23, 0.23, 1.0, -1.23, 0.23, 1.0, 1.0, -1.23, 0.23, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -0.95, -0.05, 1.0, 1.0, -0.95, -0.05, 1.0, 1.0, -0.95, -0.05, 1.0, 1.0, -0.79, -0.21, 1.0, 1.0, -0.79, -0.21, 1.0, 1.0, -0.79, -0.21, 1.0, 1.0, -0.42, -0.58, 1.0, 1.0, -0.42, -0.58, 1.0, 1.0, -0.42, -0.58, 1.0, 1.0, -0.05, -0.95, 1.0, 1.0, -0.05, -0.95, 1.0, 1.0, -0.05, -0.95, 1.0, 1.0, -0.26, -0.74, 1.0, 1.0, -0.26, -0.74, 1.0, 1.0, -0.26, -0.74, 1.0, 1.0, -3.2, 1.0, 1.0, -0.042, -0.063, 1.0, 1.0, 1.0, 1.0, 0.26, -1.26, 1.0, 1.0, -0.8, 1.0, 1.0, 0.26, -1.26, 1.0, 1.0, 0.26, -1.26, 1.0, 1.0, 0.16, -1.16, 1.0, 1.0, 0.16, -1.16, 1.0, 1.0, 0.16, -1.16, 1.0, 1.0, -0.16, -0.84, 1.0, 1.0, -0.16, -0.84, 1.0, 1.0, -0.16, -0.84, 1.0, 1.0, 1.0, 0.494, 2.27424, 0.79, -37.0, 0.506, 1.0, 1.0, 0.492, 2.2632, 0.53, -47.0, 0.508, 1.0, 0.62, -1.62, 1.0, 1.0, 0.62, -1.62, 1.0, 1.0, 0.62, -1.62, 1.0, 1.0, 0.4703, -0.0928, 1.40015, -0.247, 0.1726, -0.3122, 1.783, 1.0, 0.4273, -0.0361, 1.20404, -0.157, -0.2399, 1.0, 1.0, 0.4663, -0.0361, 1.43498, -0.157, -0.2789, 1.0, 1.0, 0.55, -0.52, 0.6, 1.0, 1.0, -1.0, 1.8, 1.0, 0.017, -0.33, 1.0, 1.0, 1.0, -0.33, 1.0, 0.2, 0.73, 1.0, 0.07, 1.0, 0.72, 1.0, 0.2, 0.08, 1.0, 1.0, 0.3675, 0.02536, 1.614, 0.25, -45.0, 0.875, 0.6325, 1.0, 1.0, 0.365, 0.02538, 1.59, 0.2, -55.0, 0.875, 0.635, 1.0, -0.828, 1.0, 0.012, -1.42, 1.0, -0.095, -0.02, -0.0467, 1.0, -0.808, 1.0, 0.0205, -1.84, -0.0022, -0.0192, 1.0, 0.679, -0.095, -0.02, -0.0467, 1.0, -0.808, 1.0, 0.0205, -1.84, -0.0022, -0.0192, 1.0, 0.679, -0.095, -0.02, -0.0467, 1.0, -1.0, -5.2, 1.0, 1.0, -1.0, -6.7, 1.0, 1.0, 1.0, 1.0, 1.0, -0.8, 1.0, 1.0, 0.482, 1.0, 2.217, 0.498, 0.02, 0.79, 1.0, 0.474, 1.0, 2.18, 0.02, 0.53, 0.506, 1.0, 0.07, 0.1, 1.0, 0.83, 1.0, 0.07, 1.0, 0.33, 0.6, 1.0, 0.33, 0.07, 1.0, 0.6, 1.0, -0.125, 0.01812, -0.65, 0.625, 1.125, 1.0, 1.0, 1.0, 1.0, 1.0, -0.25, 1.0, 0.03625, -1.36562, 1.25, 1.03125, 0.21875, -30.0, 1.0, -0.25, 1.0, 0.03625, -1.38375, 1.25, 0.21875, 1.03125, -35.0, 1.0, -0.706, 1.0, 0.0129, -1.61, 1.072, -0.027, -0.128, -0.1203, 1.0, -0.69, 1.0, 0.0195, -1.84, -0.0012, -0.0159, 1.072, 0.534, -0.027, -0.128, -0.1203, 1.0, -0.69, 1.0, 0.0195, -1.84, -0.0012, -0.0159, 1.0, 0.534, -0.027, -0.128, -0.1203, 1.0, 0.181, 1.0, 0.11, 0.709, 1.0, 0.051, 1.0, 0.055, 0.894, 1.0, 0.036, 1.0, 0.964, 1.0, -1.0, -5.3, 1.0, 1.0, -10.1, 1.0, 0.92, -0.042, -0.063, 0.08, 1.0, 1.0, 0.4, -45.0, 1.0, -1.0, -4.35, 1.0, -1.0, -2.1, 1.0, 1.0, 0.8, 1.0, 1.0, -24.0, 1.0, -64.3, 1.0, 1.0, -27.4, 1.0, 1.0, 9.1, 1.0, -0.042, -0.063, 1.0, 1.0, -0.026, -0.182, -0.1742, -0.36, -0.134, 0.826, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; 
double rhs[] = {52.6, -524.9, -0.0, -0.0, -0.0, 43.0, 44.9, -0.0, -0.0, -0.0, -0.0, 2.5, -1231.6, -0.0, 22.7, 23.4, -0.0, -0.0, 108.0, 50.0, 13.0, -2366.0, 200.0, 265.0, 300.0, -0.0, 31.0, 60.0, 134.0, 34.0, 413.0, 41.5, 15.0, 20.6, 440.0, 290.0, 3.1, 9.1, -0.0, -0.0, 34.4, 15.6, 19.2, 6.1, 13.2, -0.0, -0.0, 31.2, -1080.0, -0.0, -0.0, -0.0, 107.0, 23.4244, 23.4244, 5.288891304347825, 44.9459, 44.9459, 9.1101, 15.6166, 54.42459683471149, 385.65456097560974, 107.108, 107.108, 108.109, 108.109, 108.109, 50.051, 50.051, 50.051, 13.014, 13.014, 13.014, 200.201, 200.201, 200.201, 265.266, 265.266, 265.266, 300.301, 300.301, 300.301, 9.1101, 15.6166, 134.135, 54.42459683471149, 134.135, 134.135, 31.032, 31.032, 31.032, 60.061, 60.061, 60.061, 52.653600000000004, 52.653600000000004, 34.035, 34.035, 34.035, 23.299653954010097, 41.5425, 41.5425, 16.0, 43.539877467769195, 15.016, 15.016, 3.1041000000000003, 3.1041000000000003, 23.567400000000003, 23.567400000000003, 413.414, 413.414, 427.5267731958763, 119.16441845360825, 57.804824046007084, 31.2322, 54.42459683471149, 43.044, 43.044, 19.2202, 19.2202, 19.2202, 6.1071, 6.1071, 16.49748, 16.49748, 385.64745522388057, 385.64745522388057, 440.441, 13.2142, 13.2142, 13.2142, 155.10281253731344, 1.9692706521739132, 52.653600000000004, 452.1243262183199, 936.5421757379482, 43.044, 44.9459, 1.453505402173913, 2.1797581032608697, 9.1101, 13.5, 15.6166, 31.2322, 54.42459683471149}; 
double obj[] = {-3280.0, -3280.0, 3310.0, -1890.0, 0.0, 0.0, -1890.0, -903.0, 0.0, 432.0, 432.0, 432.0, 446.0, 446.0, 446.0, 450.0, 450.0, 450.0, 459.0, 459.0, 459.0, 483.0, 483.0, 483.0, 500.0, 500.0, 500.0, 493.0, 493.0, 493.0, -1890.0, -903.0, 506.0, 0.0, 506.0, 506.0, 505.0, 505.0, 505.0, 499.0, 499.0, 499.0, 0.0, 0.0, 512.0, 512.0, 512.0, 70.9, 39.8, 39.8, 2.04, 0.0, 1.8, 1.8, -2600.0, -2600.0, 10.4, 10.4, 28.8, 43.4, 30.4, 0.0, 0.0, -1218.0, 0.0, 0.0, 0.0, -1322.0, -1322.0, -1322.0, -1660.0, -1670.0, 14.8, 14.8, 28.8, 43.0, 30.0, -1763.0, -1722.0, -1680.0, 0.0, -1890.0, 1780.0, 1600.0, 903.0, 1760.0, 2100.0, 1000.0, 1000.0, -1890.0, 92.1, -903.0, -1218.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
