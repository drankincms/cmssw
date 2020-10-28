#include "HeterogeneousCore/SonicTriton/interface/TritonConverterBase.h"

class Int64StandardConverter : public TritonConverterBase<int64_t> {
public:
  Int64StandardConverter(const edm::ParameterSet& conf)
      : TritonConverterBase<int64_t>(conf) {}

  uint8_t* convert(const int64_t* in) const override;
};

DEFINE_EDM_PLUGIN(TritonConverterFactory<int64_t>, Int64StandardConverter, "Int64StandardConverter");

uint8_t* Int64StandardConverter::convert(const int64_t *in) {
    return reinterpret_cast<const uint8_t*>(in);
}
