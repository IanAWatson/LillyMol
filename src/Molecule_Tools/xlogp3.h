#ifndef MOLECULE_TOOLS_XLOGP3_H_
#define MOLECULE_TOOLS_XLOGP3_H_

namespace xlogp3 {

enum XlogPType {
  kNone = 0,
  kLipophilicCH3 = 1,
  kCh3XPi = 2,
  kCh3X = 3,
  kCh3Pi = 4,
  kC33h = 5,
  kC32h_lipo = 6,
  kC32hxpi = 7,
  kC32hx = 8,
  kC32hpi =  9,
  kC32h = 10,
  kC3hxpi = 11,
  kC3hx = 12,
  kC3hpi = 13,
  kC3h = 14,
  kC3xpi = 15,
  kC3x = 16,
  kC3pi = 17,
  kC3 = 18,

  kCarhx = 19,
  kCarh = 20,
  kCarar = 21,
  kCarxy = 22,
  kCary = 23,
  kCarx = 24,
  kCar = 25,

  kC22h = 26,
  kC2hcx = 27,
  kC2hcring = 28,
  kC2hc = 29,
  kC2hx = 30,
  kC2cx = 31,
  kC2cring = 32,
  kC2c = 33,
  kC2xx = 34,
  kC2xring = 35,
  kC2dx = 36,

  kC1d = 37,
  kC1 = 38,


  kNam2h = 39,
  kN32hpi = 40,
  kN32h = 41,
  kNamh = 42,
  kN3hpi = 43,
  kN3h = 44,
  kNam = 45,
  kN3pi = 46,
  kN3 = 47,
  kNarx2 = 48,
  kNarx = 49,
  kNar = 50,
  kNarhx = 51,
  kNarh = 52,
  kNarx2_2 = 53,
  kNarx_2 = 54,
  kNar_2 = 55,
  kN2h = 56,
  kN2cring = 57,
  kN2c = 58,
  kn2x = 59,

  kO3hpi = 60,
  kO3h = 61,
  kOar = 62,
  kO3pi = 63,
  kO3 = 64,
  kO2c = 65,
  kO2x = 66,
  kS3h = 67,
  kSar = 68,
  kS3x = 69,
  kS3 = 70,
  kS2c = 71,
  kS2x = 72,
  kSo = 73,
  kSo2 = 64,
  kP3 = 75,

  kFpi = 76,
  kF = 77,
  kClpi = 78,
  kCl = 79,
  kBrpi = 80,
  kBr = 81,
  kIpi = 82,
  kI = 83,

  kCyano = 84,
  kDiazo =  85,
  kNitro = 86,
  kNOxide = 87,

  kFailed = 99


};

std::optional<double> XLogP(Molecule& m);
std::optional<double> XLogP(Molecule& m, int* status);

}  // namespace xlogp3

#endif // MOLECULE_TOOLS_XLOGP3_H_
