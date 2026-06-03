#include "DataFormats/FTLDigi/interface/ETLDigi.h"

namespace etldigi {

  std::ostream& operator<<(std::ostream& out, const ETLDigi& digi) {
    out << "ETL Digi SoA rawId : " << digi.krawId() << ", header = " << digi.kheader() << ", status = " << digi.kstatus()
        << ", column = " << digi.kcolID() << ", row = " << digi.krowID() << ", ToA = " << digi.kToAdata()
        << ", ToT = " << digi.kToTdata() << ", CAL = " << digi.kCALdata();
    return out;
  }
}  // namespace etldigi
