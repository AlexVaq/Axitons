#ifndef	enumVarsGuard
	#define	enumVarsGuard

	#include <hdf5.h>

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
			FieldBase	= 1,
			FieldDev	= 2,
			FieldMisc	= 4,
			FieldBaseDev	= 3,
			FieldBaseMisc	= 5,
			FieldDevMisc	= 6,
			FieldAll	= 7,
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

		typedef	enum	FieldStatus_s {
			FieldUndefined	= 0,
			FieldGpu	= 1,
			FieldCpu	= 2,
			FieldCoherent	= 3,
		}	FieldStatus;

#ifdef	__NVCC__
	#define	Attr	inline constexpr __host__ __device__
#else
	#define	Attr	inline constexpr
#endif
		template<typename enumFlag>
		Attr enumFlag  operator &  (enumFlag  lhs, const enumFlag rhs) { return static_cast<enumFlag>(static_cast<int>(lhs) & static_cast<int>(rhs)); }
		template<typename enumFlag>
		Attr enumFlag& operator &= (enumFlag &lhs, const enumFlag rhs) { lhs  = static_cast<enumFlag>(static_cast<int>(lhs) & static_cast<int>(rhs)); return lhs; }
		template<typename enumFlag>
		Attr enumFlag  operator |  (enumFlag  lhs, const enumFlag rhs) { return static_cast<enumFlag>(static_cast<int>(lhs) | static_cast<int>(rhs)); }
		template<typename enumFlag>
		Attr enumFlag& operator |= (enumFlag &lhs, const enumFlag rhs) { lhs  = static_cast<enumFlag>(static_cast<int>(lhs) | static_cast<int>(rhs)); return lhs; }
		template<typename enumFlag>
		Attr enumFlag  operator ^  (enumFlag  lhs, const enumFlag rhs) { return static_cast<enumFlag>(static_cast<int>(lhs) ^ static_cast<int>(rhs)); }
		template<typename enumFlag>
		Attr enumFlag& operator ^= (enumFlag &lhs, const enumFlag rhs) { lhs  = static_cast<enumFlag>(static_cast<int>(lhs) ^ static_cast<int>(rhs)); return lhs; }
#undef	Attr
	}

	using namespace AxitonEnum;
#endif
