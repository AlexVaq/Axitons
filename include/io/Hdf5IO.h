#ifndef	Hdf5ReadWriterGuard
	#define	Hdf5ReadWriterGuard

	#include <string>
	#include "enum-vars.h"
	#include "cosmos/cosmos.h"
	#include "fields/fields.h"

	class	Hdf5ReadWriter {

		private:

		hid_t		file_id;
		int		fIndex;
		std::string	outputDir;
		std::string	outputName;

		public:

			 Hdf5ReadWriter	(iParms myParms);
			~Hdf5ReadWriter	();

		herr_t	readAttribute	(std::string attName, hid_t h5Type, void *data);
		herr_t	readAttribute	(std::string attName, hid_t h5Type, void *data, hid_t gid);
		herr_t	writeAttribute	(std::string attName, hid_t h5Type, void *data);
		herr_t	writeAttribute	(std::string attName, hid_t h5Type, void *data, hid_t gid);

		herr_t	readConf	(Cosmos *bck, Axiton *field);
		herr_t	writeConf	(Cosmos *bck, Axiton *field);

		herr_t	writeData	(std::string dataPath, hid_t h5Type, void *data, size_t aSize);
		herr_t	writeData	(std::string dataPath, hid_t h5Type, void *data, size_t aSize, hid_t gid);

		herr_t	nextFile	(int jump = 0);
		hid_t	currentFile	() { return file_id; }
		hid_t	currentIndex	() { return fIndex;  }
	};
#endif
