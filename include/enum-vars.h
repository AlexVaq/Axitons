#ifndef	enumVarsGuard
	#define	enumVarsGuard

	typedef	unsigned int uint;

	namespace AxitonEnum {
		typedef enum	FieldPrecision_s {
			SinglePrecision = 4,
			DoublePrecision = 8,
		}	FieldPrecision;

		typedef enum	AllocType_s {
			TrackHostAlloc,
			TrackDeviceAlloc,
		}	AllocType;

		typedef	enum	MemDirection_s {
			DeviceToHost,
			HostToDevice,
		}	MemDirection;

		typedef	enum	FieldIndex_s {
			FieldBase,
			FieldDev,
			FieldMisc,
		}	FieldIndex;

		typedef	enum	FieldExpansion_s {
			Minkowski,
			Radiation,
		}	FieldExpansion;

		typedef	enum	FieldType_s {
			FieldCompact,
			FieldNonCompact,
		}	FieldType;

		typedef enum	InitialCond_s {
			IcNone,
			IcRead,
			IcFlat,
			IcSinc2,
		}	InitialCond;

		typedef enum	PropagatorType_s {
			PropagatorLeapFrog2,
			PropagatorLeapFrog4,
			PropagatorOmelyan2,
			PropagatorOmelyan4,
			PropagatorRKN4,
		}	PropagatorType;
	}

	using namespace AxitonEnum;
#endif
