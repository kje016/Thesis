from sage.all import *

BG2= [
{0: [9, 174, 0, 72, 3, 156, 143, 145], 1: [117, 97, 0, 110, 26, 143, 19, 131], 2: [204, 166, 0, 23, 53, 14, 176, 71], 3: [26, 66, 0, 181, 35, 3, 165, 21],6: [189, 71, 0, 95, 115, 40, 196, 23], 9: [205, 172, 0, 8, 127, 123, 13, 112],10: [0, 0, 0, 1, 0, 0, 0, 1], 11: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [167, 27, 137, 53, 19, 17, 18, 142], 3: [166, 36, 124, 156, 94, 65, 27, 174], 4: [253, 48, 0, 115, 104, 63, 3, 183], 5: [125, 92, 0, 156, 66, 1, 102, 27], 6: [226, 31, 88, 115, 84, 55, 185, 96], 7: [156, 187, 0, 200, 98, 37, 17, 23], 8: [224, 185, 0, 29, 69, 171, 14, 9], 9: [252, 3, 55, 31, 50, 133, 180, 167], 11: [0, 0, 0, 0, 0, 0, 0, 0], 12: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [81, 25, 20, 152, 95, 98, 126, 74], 1: [114, 114, 94, 131, 106, 168, 163, 31], 3: [44, 117, 99, 46, 92, 107, 47, 3], 4: [52, 110, 9, 191, 110, 82, 183, 53], 8: [240, 114, 108, 91, 111, 142, 132, 155], 10: [1, 1, 1, 0, 1, 1, 1, 0], 12: [0, 0, 0, 0, 0, 0, 0, 0], 13: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [8, 136, 38, 185, 120, 53, 36, 239], 2: [58, 175, 15, 6, 121, 174, 48, 171], 4: [158, 113, 102, 36, 22, 174, 18, 95], 5: [104, 72, 146, 124, 4, 127, 111, 110], 6: [209, 123, 12, 124, 73, 17, 203, 159], 7: [54, 118, 57, 110, 49, 89, 3, 199], 8: [18, 28, 53, 156, 128, 17, 191, 43], 9: [128, 186, 46, 133, 79, 105, 160, 75], 10: [0, 0, 0, 1, 0, 0, 0, 1], 13: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [179, 72, 0, 200, 42, 86, 43, 29], 1: [214, 74, 136, 16, 24, 67, 27, 140], 11: [71, 29, 157, 101, 51, 83, 117, 180], 14: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [231, 10, 0, 185, 40, 79, 136, 121], 1: [41, 44, 131, 138, 140, 84, 49, 41], 5: [194, 121, 142, 170, 84, 35, 36, 169], 7: [159, 80, 141, 219, 137, 103, 132, 88], 11: [103, 48, 64, 193, 71, 60, 62, 207], 15: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [155, 129, 0, 123, 109, 47, 7, 137], 5: [228, 92, 124, 55, 87, 154, 34, 72], 7: [45, 100, 99, 31, 107, 10, 198, 172], 9: [28, 49, 45, 222, 133, 155, 168, 124], 11: [158, 184, 148, 209, 139, 29, 12, 56], 16: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [129, 80, 0, 103, 97, 48, 163, 86], 5: [147, 186, 45, 13, 135, 125, 78, 186], 7: [140, 16, 148, 105, 35, 24, 143, 87], 11: [3, 102, 96, 150, 108, 47, 107, 172], 13: [116, 143, 78, 181, 65, 55, 58, 154], 17: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [142, 118, 0, 147, 70, 53, 101, 176], 1: [94, 70, 65, 43, 69, 31, 177, 169], 12: [230, 152, 87, 152, 88, 161, 22, 225], 18: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [203, 28, 0, 2, 97, 104, 186, 167], 8: [205, 132, 97, 30, 40, 142, 27, 238], 10: [61, 185, 51, 184, 24, 99, 205, 48], 11: [247, 178, 85, 83, 49, 64, 81, 68], 19: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [11, 59, 0, 174, 46, 111, 125, 38], 1: [185, 104, 17, 150, 41, 25, 60, 217], 6: [0, 22, 156, 8, 101, 174, 177, 208], 7: [117, 52, 20, 56, 96, 23, 51, 232], 20: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [11, 32, 0, 99, 28, 91, 39, 178], 7: [236, 92, 7, 138, 30, 175, 29, 214], 9: [210, 174, 4, 110, 116, 24, 35, 168], 13: [56, 154, 2, 99, 64, 141, 8, 51], 21: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [63, 39, 0, 46, 33, 122, 18, 124], 3: [111, 93, 113, 217, 122, 11, 155, 122], 11: [14, 11, 48, 109, 131, 4, 49, 72], 22: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [83, 49, 0, 37, 76, 29, 32, 48], 1: [2, 125, 112, 113, 37, 91, 53, 57], 8: [38, 35, 102, 143, 62, 27, 95, 167], 13: [222, 166, 26, 140, 47, 127, 186, 219], 23: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [115, 19, 0, 36, 143, 11, 91, 82], 6: [145, 118, 138, 95, 51, 145, 20, 232], 11: [3, 21, 57, 40, 130, 8, 52, 204], 13: [232, 163, 27, 116, 97, 166, 109, 162], 24: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [51, 68, 0, 116, 139, 137, 174, 38], 10: [175, 63, 73, 200, 96, 103, 108, 217], 11: [213, 81, 99, 110, 128, 40, 102, 157], 25: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [203, 87, 0, 75, 48, 78, 125, 170], 9: [142, 177, 79, 158, 9, 158, 31, 23], 11: [8, 135, 111, 134, 28, 17, 54, 175], 12: [242, 64, 143, 97, 8, 165, 176, 202], 26: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [254, 158, 0, 48, 120, 134, 57, 196], 5: [124, 23, 24, 132, 43, 23, 201, 173], 11: [114, 9, 109, 206, 65, 62, 142, 195], 12: [64, 6, 18, 2, 42, 163, 35, 218], 27: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [220, 186, 0, 68, 17, 173, 129, 128], 6: [194, 6, 18, 16, 106, 31, 203, 211], 7: [50, 46, 86, 156, 142, 22, 140, 210], 28: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [87, 58, 0, 35, 79, 13, 110, 39], 1: [20, 42, 158, 138, 28, 135, 124, 84], 10: [185, 156, 154, 86, 41, 145, 52, 88], 29: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [26, 76, 0, 6, 2, 128, 196, 117], 4: [105, 61, 148, 20, 103, 52, 35, 227], 11: [29, 153, 104, 141, 78, 173, 114, 6], 30: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [76, 157, 0, 80, 91, 156, 10, 238], 8: [42, 175, 17, 43, 75, 166, 122, 13], 13: [210, 67, 33, 81, 81, 40, 23, 11], 31: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [222, 20, 0, 49, 54, 18, 202, 195], 2: [63, 52, 4, 1, 132, 163, 126, 44], 32: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [23, 106, 0, 156, 68, 110, 52, 5], 3: [235, 86, 75, 54, 115, 132, 170, 94], 5: [238, 95, 158, 134, 56, 150, 13, 111], 33: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [46, 182, 0, 153, 30, 113, 113, 81], 2: [139, 153, 69, 88, 42, 108, 161, 19], 9: [8, 64, 87, 63, 101, 61, 88, 130], 34: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [228, 45, 0, 211, 128, 72, 197, 66], 5: [156, 21, 65, 94, 63, 136, 194, 95], 35: [0, 0, 0, 0, 0, 0, 0, 0]},
{2: [29, 67, 0, 90, 142, 36, 164, 146], 7: [143, 137, 100, 6, 28, 38, 172, 66], 12: [160, 55, 13, 221, 100, 53, 49, 190], 13: [122, 85, 7, 6, 133, 145, 161, 86], 36: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [8, 103, 0, 27, 13, 42, 168, 64], 6: [151, 50, 32, 118, 10, 104, 193, 181], 37: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [98, 70, 0, 216, 106, 64, 14, 7], 2: [101, 111, 126, 212, 77, 24, 186, 144], 5: [135, 168, 110, 193, 43, 149, 46, 16], 38: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [18, 110, 0, 108, 133, 139, 50, 25], 4: [28, 17, 154, 61, 25, 161, 27, 57], 39: [0, 0, 0, 0, 0, 0, 0, 0]},
{2: [71, 120, 0, 106, 87, 84, 70, 37], 5: [240, 154, 35, 44, 56, 173, 17, 139], 7: [9, 52, 51, 185, 104, 93, 50, 221], 9: [84, 56, 134, 176, 70, 29, 6, 17], 40: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [106, 3, 0, 147, 80, 117, 115, 201], 13: [1, 170, 20, 182, 139, 148, 189, 46], 41: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [242, 84, 0, 108, 32, 116, 110, 179], 5: [44, 8, 20, 21, 89, 73, 0, 14], 12: [166, 17, 122, 110, 71, 142, 163, 116], 42: [0, 0, 0, 0, 0, 0, 0, 0]},
{2: [132, 165, 0, 71, 135, 105, 163, 46], 7: [164, 179, 88, 12, 6, 137, 173, 2], 10: [235, 124, 13, 109, 2, 29, 179, 106], 43: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [147, 173, 0, 29, 37, 11, 197, 184], 12: [85, 177, 19, 201, 25, 41, 191, 135], 13: [36, 12, 78, 69, 114, 162, 193, 141], 44: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [57, 77, 0, 91, 60, 126, 157, 85], 5: [40, 184, 157, 165, 137, 152, 167, 225], 11: [63, 18, 6, 55, 93, 172, 181, 175], 45: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [140, 25, 0, 1, 121, 73, 197, 178], 2: [38, 151, 63, 175, 129, 154, 167, 112], 7: [154, 170, 82, 83, 26, 129, 179, 106], 46: [0, 0, 0, 0, 0, 0, 0, 0]},
{10: [219, 37, 0, 40, 97, 167, 181, 154], 13: [151, 31, 144, 12, 56, 38, 193, 114], 47: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [31, 84, 0, 37, 1, 112, 157, 42], 5: [66, 151, 93, 97, 70, 7, 173, 41], 11: [38, 190, 19, 46, 1, 19, 191, 105], 48: [0, 0, 0, 0, 0, 0, 0, 0]},
{0: [239, 93, 0, 106, 119, 109, 181, 167], 7: [172, 132, 24, 181, 32, 6, 157, 45], 12: [34, 57, 138, 154, 142, 105, 173, 189], 49: [0, 0, 0, 0, 0, 0, 0, 0]},
{2: [0, 103, 0, 98, 6, 160, 193, 78], 10: [75, 107, 36, 35, 73, 156, 163, 67], 13: [120, 163, 143, 36, 102, 82, 179, 180], 50: [0, 0, 0, 0, 0, 0, 0, 0]},
{1: [129, 147, 0, 120, 48, 132, 191, 53], 5: [229, 7, 2, 101, 47, 6, 197, 215], 11: [118, 60, 55, 81, 19, 8, 167, 230], 51: [0, 0, 0, 0, 0, 0, 0, 0]}
]

BG1 = [{0: [250, 307, 73, 223, 211, 294, 0, 135], 1: [69, 19, 15, 16, 198, 118, 0, 227], 2: [226, 50, 103, 94, 188, 167, 0, 126], 3: [159, 369, 49, 91, 186, 330, 0, 134], 5: [100, 181, 240, 74, 219, 207, 0, 84], 6: [10, 216, 39, 10, 4, 165, 0, 83], 9: [59, 317, 15, 0, 29, 243, 0, 53], 10: [229, 288, 162, 205, 144, 250, 0, 225], 11: [110, 109, 215, 216, 116, 1, 0, 205], 12: [191, 17, 164, 21, 216, 339, 0, 128], 13: [9, 357, 133, 215, 115, 201, 0, 75], 15: [195, 215, 298, 14, 233, 53, 0, 135], 16: [23, 106, 110, 70, 144, 347, 0, 217], 18: [190, 242, 113, 141, 95, 304, 0, 220], 19: [35, 180, 16, 198, 216, 167, 0, 90], 20: [239, 330, 189, 104, 73, 47, 0, 105], 21: [31, 346, 32, 81, 261, 188, 0, 137], 22: [1, 1, 1, 1, 1, 1, 0, 1], 23: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [2, 76, 303, 141, 179, 77, 22, 96], 2: [239, 76, 294, 45, 162, 225, 11, 236], 3: [117, 73, 27, 151, 223, 96, 124, 136], 4: [124, 288, 261, 46, 256, 338, 0, 221], 5: [71, 144, 161, 119, 160, 268, 10, 128], 7: [222, 331, 133, 157, 76, 112, 0, 92], 8: [104, 331, 4, 133, 202, 302, 0, 172], 9: [173, 178, 80, 87, 117, 50, 2, 56], 11: [220, 295, 129, 206, 109, 167, 16, 11], 12: [102, 342, 300, 93, 15, 253, 60, 189], 14: [109, 217, 76, 79, 72, 334, 0, 95], 15: [132, 99, 266, 9, 152, 242, 6, 85], 16: [142, 354, 72, 118, 158, 257, 30, 153], 17: [155, 114, 83, 194, 147, 133, 0, 87], 19: [255, 331, 260, 31, 156, 9, 168, 163], 21: [28, 112, 301, 187, 119, 302, 31, 216], 22: [0, 0, 0, 0, 0, 0, 105, 0], 23: [0, 0, 0, 0, 0, 0, 0, 0], 24: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [106, 205, 68, 207, 258, 226, 132, 189], 1: [111, 250, 7, 203, 167, 35, 37, 4], 2: [185, 328, 80, 31, 220, 213, 21, 225], 4: [63, 332, 280, 176, 133, 302, 180, 151], 5: [117, 256, 38, 180, 243, 111, 4, 236], 6: [93, 161, 227, 186, 202, 265, 149, 117], 7: [229, 267, 202, 95, 218, 128, 48, 179], 8: [177, 160, 200, 153, 63, 237, 38, 92], 9: [95, 63, 71, 177, 0, 294, 122, 24], 10: [39, 129, 106, 70, 3, 127, 195, 68], 13: [142, 200, 295, 77, 74, 110, 155, 6], 14: [225, 88, 283, 214, 229, 286, 28, 101], 15: [225, 53, 301, 77, 0, 125, 85, 33], 17: [245, 131, 184, 198, 216, 131, 47, 96], 18: [205, 240, 246, 117, 269, 163, 179, 125], 19: [251, 205, 230, 223, 200, 210, 42, 67], 20: [117, 13, 276, 90, 234, 7, 66, 230], 24: [0, 0, 0, 0, 0, 0, 0, 0], 25: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [121, 276, 220, 201, 187, 97, 4, 128], 1: [89, 87, 208, 18, 145, 94, 6, 23], 3: [84, 0, 30, 165, 166, 49, 33, 162], 4: [20, 275, 197, 5, 108, 279, 113, 220], 6: [150, 199, 61, 45, 82, 139, 49, 43], 7: [131, 153, 175, 142, 132, 166, 21, 186], 8: [243, 56, 79, 16, 197, 91, 6, 96], 10: [136, 132, 281, 34, 41, 106, 151, 1], 11: [86, 305, 303, 155, 162, 246, 83, 216], 12: [246, 231, 253, 213, 57, 345, 154, 22], 13: [219, 341, 164, 147, 36, 269, 87, 24], 14: [211, 212, 53, 69, 115, 185, 5, 167], 16: [240, 304, 44, 96, 242, 249, 92, 200], 17: [76, 300, 28, 74, 165, 215, 173, 32], 18: [244, 271, 77, 99, 0, 143, 120, 235], 20: [144, 39, 319, 30, 113, 121, 2, 172], 21: [12, 357, 68, 158, 108, 121, 142, 219], 22: [1, 1, 1, 1, 1, 1, 0, 1], 25: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [157, 332, 233, 170, 246, 42, 24, 64], 1: [102, 181, 205, 10, 235, 256, 204, 211], 26: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [205, 195, 83, 164, 261, 219, 185, 2], 1: [236, 14, 292, 59, 181, 130, 100, 171], 3: [194, 115, 50, 86, 72, 251, 24, 47], 12: [231, 166, 318, 80, 283, 322, 65, 143], 16: [28, 241, 201, 182, 254, 295, 207, 210], 21: [123, 51, 267, 130, 79, 258, 161, 180], 22: [115, 157, 279, 153, 144, 283, 72, 180], 27: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [183, 278, 289, 158, 80, 294, 6, 199], 6: [22, 257, 21, 119, 144, 73, 27, 22], 10: [28, 1, 293, 113, 169, 330, 163, 23], 11: [67, 351, 13, 21, 90, 99, 50, 100], 13: [244, 92, 232, 63, 59, 172, 48, 92], 17: [11, 253, 302, 51, 177, 150, 24, 207], 18: [157, 18, 138, 136, 151, 284, 38, 52], 20: [211, 225, 235, 116, 108, 305, 91, 13], 28: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [220, 9, 12, 17, 169, 3, 145, 77], 1: [44, 62, 88, 76, 189, 103, 88, 146], 4: [159, 316, 207, 104, 154, 224, 112, 209], 7: [31, 333, 50, 100, 184, 297, 153, 32], 8: [167, 290, 25, 150, 104, 215, 159, 166], 14: [104, 114, 76, 158, 164, 39, 76, 18], 29: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [112, 307, 295, 33, 54, 348, 172, 181], 1: [4, 179, 133, 95, 0, 75, 2, 105], 3: [7, 165, 130, 4, 252, 22, 131, 141], 12: [211, 18, 231, 217, 41, 312, 141, 223], 16: [102, 39, 296, 204, 98, 224, 96, 177], 19: [164, 224, 110, 39, 46, 17, 99, 145], 21: [109, 368, 269, 58, 15, 59, 101, 199], 22: [241, 67, 245, 44, 230, 314, 35, 153], 24: [90, 170, 154, 201, 54, 244, 116, 38], 30: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [103, 366, 189, 9, 162, 156, 6, 169], 1: [182, 232, 244, 37, 159, 88, 10, 12], 10: [109, 321, 36, 213, 93, 293, 145, 206], 11: [21, 133, 286, 105, 134, 111, 53, 221], 13: [142, 57, 151, 89, 45, 92, 201, 17], 17: [14, 303, 267, 185, 132, 152, 4, 212], 18: [61, 63, 135, 109, 76, 23, 164, 92], 20: [216, 82, 209, 218, 209, 337, 173, 205], 31: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [98, 101, 14, 82, 178, 175, 126, 116], 2: [149, 339, 80, 165, 1, 253, 77, 151], 4: [167, 274, 211, 174, 28, 27, 156, 70], 7: [160, 111, 75, 19, 267, 231, 16, 230], 8: [49, 383, 161, 194, 234, 49, 12, 115], 14: [58, 354, 311, 103, 201, 267, 70, 84], 32: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [77, 48, 16, 52, 55, 25, 184, 45], 1: [41, 102, 147, 11, 23, 322, 194, 115], 12: [83, 8, 290, 2, 274, 200, 123, 134], 16: [182, 47, 289, 35, 181, 351, 16, 1], 21: [78, 188, 177, 32, 273, 166, 104, 152], 22: [252, 334, 43, 84, 39, 338, 109, 165], 23: [22, 115, 280, 201, 26, 192, 124, 107], 33: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [160, 77, 229, 142, 225, 123, 6, 186], 1: [42, 186, 235, 175, 162, 217, 20, 215], 10: [21, 174, 169, 136, 244, 142, 203, 124], 11: [32, 232, 48, 3, 151, 110, 153, 180], 13: [234, 50, 105, 28, 238, 176, 104, 98], 18: [7, 74, 52, 182, 243, 76, 207, 80], 34: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [177, 313, 39, 81, 231, 311, 52, 220], 3: [248, 177, 302, 56, 0, 251, 147, 185], 7: [151, 266, 303, 72, 216, 265, 1, 154], 20: [185, 115, 160, 217, 47, 94, 16, 178], 23: [62, 370, 37, 78, 36, 81, 46, 150], 35: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [206, 142, 78, 14, 0, 22, 1, 124], 12: [55, 248, 299, 175, 186, 322, 202, 144], 15: [206, 137, 54, 211, 253, 277, 118, 182], 16: [127, 89, 61, 191, 16, 156, 130, 95], 17: [16, 347, 179, 51, 0, 66, 1, 72], 21: [229, 12, 258, 43, 79, 78, 2, 76], 36: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [40, 241, 229, 90, 170, 176, 173, 39], 1: [96, 2, 290, 120, 0, 348, 6, 138], 10: [65, 210, 60, 131, 183, 15, 81, 220], 13: [63, 318, 130, 209, 108, 81, 182, 173], 18: [75, 55, 184, 209, 68, 176, 53, 142], 25: [179, 269, 51, 81, 64, 113, 46, 49], 37: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [64, 13, 69, 154, 270, 190, 88, 78], 3: [49, 338, 140, 164, 13, 293, 198, 152], 11: [49, 57, 45, 43, 99, 332, 160, 84], 20: [51, 289, 115, 189, 54, 331, 122, 5], 22: [154, 57, 300, 101, 0, 114, 182, 205], 38: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [7, 260, 257, 56, 153, 110, 91, 183], 14: [164, 303, 147, 110, 137, 228, 184, 112], 16: [59, 81, 128, 200, 0, 247, 30, 106], 17: [1, 358, 51, 63, 0, 116, 3, 219], 21: [144, 375, 228, 4, 162, 190, 155, 129], 39: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [42, 130, 260, 199, 161, 47, 1, 183], 12: [233, 163, 294, 110, 151, 286, 41, 215], 13: [8, 280, 291, 200, 0, 246, 167, 180], 18: [155, 132, 141, 143, 241, 181, 68, 143], 19: [147, 4, 295, 186, 144, 73, 148, 14], 40: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [60, 145, 64, 8, 0, 87, 12, 179], 1: [73, 213, 181, 6, 0, 110, 6, 108], 7: [72, 344, 101, 103, 118, 147, 166, 159], 8: [127, 242, 270, 198, 144, 258, 184, 138], 10: [224, 197, 41, 8, 0, 204, 191, 196], 41: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [151, 187, 301, 105, 265, 89, 6, 77], 3: [186, 206, 162, 210, 81, 65, 12, 187], 9: [217, 264, 40, 121, 90, 155, 15, 203], 11: [47, 341, 130, 214, 144, 244, 5, 167], 22: [160, 59, 10, 183, 228, 30, 30, 130], 42: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [249, 205, 79, 192, 64, 162, 6, 197], 5: [121, 102, 175, 131, 46, 264, 86, 122], 16: [109, 328, 132, 220, 266, 346, 96, 215], 20: [131, 213, 283, 50, 9, 143, 42, 65], 21: [171, 97, 103, 106, 18, 109, 199, 216], 43: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [64, 30, 177, 53, 72, 280, 44, 25], 12: [142, 11, 20, 0, 189, 157, 58, 47], 13: [188, 233, 55, 3, 72, 236, 130, 126], 17: [158, 22, 316, 148, 257, 113, 131, 178], 44: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [156, 24, 249, 88, 180, 18, 45, 185], 2: [147, 89, 50, 203, 0, 6, 18, 127], 10: [170, 61, 133, 168, 0, 181, 132, 117], 18: [152, 27, 105, 122, 165, 304, 100, 199], 45: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [112, 298, 289, 49, 236, 38, 9, 32], 3: [86, 158, 280, 157, 199, 170, 125, 178], 4: [236, 235, 110, 64, 0, 249, 191, 2], 11: [116, 339, 187, 193, 266, 288, 28, 156], 22: [222, 234, 281, 124, 0, 194, 6, 58], 46: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [23, 72, 172, 1, 205, 279, 4, 27], 6: [136, 17, 295, 166, 0, 255, 74, 141], 7: [116, 383, 96, 65, 0, 111, 16, 11], 14: [182, 312, 46, 81, 183, 54, 28, 181], 47: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [195, 71, 270, 107, 0, 325, 21, 163], 2: [243, 81, 110, 176, 0, 326, 142, 131], 4: [215, 76, 318, 212, 0, 226, 192, 169], 15: [61, 136, 67, 127, 277, 99, 197, 98], 48: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [25, 194, 210, 208, 45, 91, 98, 165], 6: [104, 194, 29, 141, 36, 326, 140, 232], 8: [194, 101, 304, 174, 72, 268, 22, 9], 49: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [128, 222, 11, 146, 275, 102, 4, 32], 4: [165, 19, 293, 153, 0, 1, 1, 43], 19: [181, 244, 50, 217, 155, 40, 40, 200], 21: [63, 274, 234, 114, 62, 167, 93, 205], 50: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [86, 252, 27, 150, 0, 273, 92, 232], 14: [236, 5, 308, 11, 180, 104, 136, 32], 18: [84, 147, 117, 53, 0, 243, 106, 118], 25: [6, 78, 29, 68, 42, 107, 6, 103], 51: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [216, 159, 91, 34, 0, 171, 2, 170], 10: [73, 229, 23, 130, 90, 16, 88, 199], 13: [120, 260, 105, 210, 252, 95, 112, 26], 24: [9, 90, 135, 123, 173, 212, 20, 105], 52: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [95, 100, 222, 175, 144, 101, 4, 73], 7: [177, 215, 308, 49, 144, 297, 49, 149], 22: [172, 258, 66, 177, 166, 279, 125, 175], 25: [61, 256, 162, 128, 19, 222, 194, 108], 53: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [221, 102, 210, 192, 0, 351, 6, 103], 12: [112, 201, 22, 209, 211, 265, 126, 110], 14: [199, 175, 271, 58, 36, 338, 63, 151], 24: [121, 287, 217, 30, 162, 83, 20, 211], 54: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [2, 323, 170, 114, 0, 56, 10, 199], 2: [187, 8, 20, 49, 0, 304, 30, 132], 11: [41, 361, 140, 161, 76, 141, 6, 172], 21: [211, 105, 33, 137, 18, 101, 92, 65], 55: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [127, 230, 187, 82, 197, 60, 4, 161], 7: [167, 148, 296, 186, 0, 320, 153, 237], 15: [164, 202, 5, 68, 108, 112, 197, 142], 17: [159, 312, 44, 150, 0, 54, 155, 180], 56: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [161, 320, 207, 192, 199, 100, 4, 231], 6: [197, 335, 158, 173, 278, 210, 45, 174], 12: [207, 2, 55, 26, 0, 195, 168, 145], 22: [103, 266, 285, 187, 205, 268, 185, 100], 57: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [37, 210, 259, 222, 216, 135, 6, 11], 14: [105, 313, 179, 157, 16, 15, 200, 207], 15: [51, 297, 178, 0, 0, 35, 177, 42], 18: [120, 21, 160, 6, 0, 188, 43, 100], 58: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [198, 269, 298, 81, 72, 319, 82, 59], 13: [220, 82, 15, 195, 144, 236, 2, 204], 23: [122, 115, 115, 138, 0, 85, 135, 161], 59: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [167, 185, 151, 123, 190, 164, 91, 121], 9: [151, 177, 179, 90, 0, 196, 64, 90], 10: [157, 289, 64, 73, 0, 209, 198, 26], 12: [163, 214, 181, 10, 0, 246, 100, 140], 60: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [173, 258, 102, 12, 153, 236, 4, 115], 3: [139, 93, 77, 77, 0, 264, 28, 188], 7: [149, 346, 192, 49, 165, 37, 109, 168], 19: [0, 297, 208, 114, 117, 272, 188, 52], 61: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [157, 175, 32, 67, 216, 304, 10, 4], 8: [137, 37, 80, 45, 144, 237, 84, 103], 17: [149, 312, 197, 96, 2, 135, 12, 30], 62: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [167, 52, 154, 23, 0, 123, 2, 53], 3: [173, 314, 47, 215, 0, 77, 75, 189], 9: [139, 139, 124, 60, 0, 25, 142, 215], 18: [151, 288, 207, 167, 183, 272, 128, 24], 63: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [149, 113, 226, 114, 27, 288, 163, 222], 4: [157, 14, 65, 91, 0, 83, 10, 170], 24: [137, 218, 126, 78, 35, 17, 162, 71], 64: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [151, 113, 228, 206, 52, 210, 1, 22], 16: [163, 132, 69, 22, 243, 3, 163, 127], 18: [173, 114, 176, 134, 0, 53, 99, 49], 25: [139, 168, 102, 161, 270, 167, 98, 125], 65: [0, 0, 0, 0, 0, 0, 0, 0]},
       {0: [139, 80, 234, 84, 18, 79, 4, 191], 7: [157, 78, 227, 4, 0, 244, 6, 211], 9: [163, 163, 259, 9, 0, 293, 142, 187], 22: [173, 274, 260, 12, 57, 272, 3, 148], 66: [0, 0, 0, 0, 0, 0, 0, 0]},
       {1: [149, 135, 101, 184, 168, 82, 181, 177], 6: [151, 149, 228, 121, 0, 67, 45, 114], 10: [167, 15, 126, 29, 144, 235, 153, 93], 67: [0, 0, 0, 0, 0, 0, 0, 0]},]

def create_BG(Zc, iLS, bg):
    M, N = 42 if bg==2 else 46, 52 if bg==2 else 68
    BG = Matrix(ZZ, M, N)
    bgm = BG1 if bg==1 else BG2
    non_getter = [-1]*Zc
    for i in range(M):
        for j in range(N):
            if bgm[i].get(j):
                BG[i,j] = bgm[i].get(j)[iLS] % Zc
            else:
                BG[i,j] = -1
    if bg == 2:
        if BG[2][N-M] == 0:
            Bi = 1
        else:
            Bi = 0
    else:
        if BG[0][22] == 1:
            Bi = (0, 1)
        else:
            Bi = (1, BG[1][22])
    return BG, Bi
