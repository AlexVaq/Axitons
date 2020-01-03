#include <cstdio>
#include "io/Hdf5IO.h"

	Hdf5ReadWriter::Hdf5ReadWriter (iParms myParms) : outputDir(myParms.outputDir), outputName(myParms.outputName), fIndex(fIndex) {

	char cIndex[8];

	sprintf (cIndex, ".%05d\0", fIndex);
	std::string	fileName(outputDir + outputName + cIndex);

        /*      Create the file and release the plist   */
	if ((file_id = H5Fcreate (fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0)
	{
		printf ("Error opening file %s", fileName.c_str());
		return;
	}
}

	Hdf5ReadWriter::~Hdf5ReadWriter () {

	H5Fclose(file_id);
}

herr_t	Hdf5ReadWriter::nextFile	(int jump=0) {

	H5Fclose(file_id);

	fIndex++;

	char cIndex[8];

	sprintf (cIndex, ".%05d\0", fIndex);
	std::string	fileName(outputDir + outputName + cIndex);

	if ((file_id = H5Fcreate (fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0)
	{
		printf ("Error opening file %s", fileName.c_str());
		return	-1;
	}

	return	0;
}

herr_t  Hdf5ReadWriter::readAttribute	(std::string attName, hid_t h5Type, void *data)
{
        hid_t   attr;
        herr_t  status;

        if ((attr   = H5Aopen_by_name (file_id, ".", attName.c_str(), H5P_DEFAULT, H5P_DEFAULT)) < 0)
                printf ("Error opening attribute %s", attName.c_str());
        if ((status = H5Aread (attr, h5Type, data)) < 0)
                printf ("Error reading attribute %s", attName.c_str());
        status = H5Aclose(attr);

        return  status;
}


herr_t  Hdf5ReadWriter::writeAttribute	(std::string attName, hid_t h5Type, void *data)
{
        hid_t   attr, attr_id;
        herr_t  status;

        attr_id = H5Screate(H5S_SCALAR);
        if ((attr   = H5Acreate2 (file_id, attName.c_str(), h5Type, attr_id, H5P_DEFAULT, H5P_DEFAULT)) < 0)
                printf ("Error creating attribute %s", attName.c_str());
        if ((status = H5Awrite (attr, h5Type, data)) < 0)
                printf ("Error writing attribute %s to file", attName.c_str());
        H5Sclose(attr_id);
        status = H5Aclose(attr);

        return  status;
}


herr_t	Hdf5ReadWriter::writeData (std::string dataName, hid_t h5Type, void *aData, size_t aSize)
{
	hid_t   base_id, dataSpace, sSpace, dataSet;
	hsize_t dims[1] = { aSize };

	/*      Create the group for the data if it doesn't exist       */
	auto status = H5Lexists (file_id, dataName.c_str(), H5P_DEFAULT);  // Create group if it doesn't exists

	if (!status)
		base_id = H5Gcreate2(file_id, dataName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	else {
		if (status > 0)
			base_id = H5Gopen2(file_id, dataName.c_str(), H5P_DEFAULT);
		else {
			printf ("Error: can't check whether group %s exists", dataName.c_str());
			return	-1;
		}
	}

	/*      Create dataset  */
	dataSpace = H5Screate_simple(1, dims, NULL);
	dataSet   = H5Dcreate(base_id, "data", h5Type, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	sSpace    = H5Dget_space (dataSet);

	/*      Write spectrum data     */
	if (H5Dwrite(dataSet, h5Type, dataSpace, sSpace, H5P_DEFAULT, aData) < 0) {
		printf ("Error writing %zu elements to dataset", aSize);
		return	-1;
	}

	/*      Close everything        */
	H5Sclose (sSpace);
	H5Dclose (dataSet);
	H5Sclose (dataSpace);
	H5Gclose (base_id);
}


herr_t	Hdf5ReadWriter::readConf  (Cosmos *bck, Axiton *field) { return	0; }
herr_t	Hdf5ReadWriter::writeConf (Cosmos *bck, Axiton *field) { return	0;}
