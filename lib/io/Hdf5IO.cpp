#include <cstdio>
#include <string>
#include "io/Hdf5IO.h"

struct	Hdf5Header {
	double		lSize;
	double		nQcd;
	double		z;
	double		R;
	std::string	exp;
	std::string	icType;
	double		icParm1;
	int		icParm2;
	std::string	pType;
	std::string	prec;
	int		nSize;
	double		wallTime;
	double massA;
};

hid_t	createGroup	(std::string gName, hid_t base_id) {

	hid_t	group_id;

	auto status = H5Lexists (base_id, gName.c_str(), H5P_DEFAULT);  // Create group if and only if it doesn't exists

	if (!status)
		group_id = H5Gcreate2(base_id, gName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	else {
		if (status > 0)
			group_id = H5Gopen2(base_id, gName.c_str(), H5P_DEFAULT);
		else {
			printf ("Error: can't check whether group %s exists\n", gName.c_str());
			return	-1;
		}
	}

	return	group_id;
}

herr_t	writeHeader	(Hdf5ReadWriter &IOHandler, Hdf5Header &myHead) {

	hid_t	gid, hdf5String;

        hdf5String = H5Tcopy(H5T_C_S1);
        H5Tset_size   (hdf5String, 32);
        H5Tset_strpad (hdf5String, H5T_STR_NULLTERM);

	gid = createGroup ("Physical",	 IOHandler.currentFile());

	IOHandler.writeAttribute("Physical Length", H5T_NATIVE_DOUBLE, &myHead.lSize,        gid);
	IOHandler.writeAttribute("nQcd",            H5T_NATIVE_DOUBLE, &myHead.nQcd,         gid);
	IOHandler.writeAttribute("z",               H5T_NATIVE_DOUBLE, &myHead.z,            gid);
	IOHandler.writeAttribute("R",               H5T_NATIVE_DOUBLE, &myHead.R,            gid);
	IOHandler.writeAttribute("Axion mass",      H5T_NATIVE_DOUBLE, &myHead.massA,        gid);
	IOHandler.writeAttribute("Expansion",       hdf5String,        &myHead.exp.at(0),    gid);

	H5Gclose (gid);

	gid = createGroup ("IC",	 IOHandler.currentFile());

	IOHandler.writeAttribute("IcType",          hdf5String,        &myHead.icType.at(0), gid);
	IOHandler.writeAttribute("Ic Parm1",        H5T_NATIVE_DOUBLE, &myHead.icParm1,      gid);
	IOHandler.writeAttribute("Ic Parm2",        H5T_NATIVE_INT,    &myHead.icParm2,      gid);

	H5Gclose (gid);

	gid = createGroup ("Propagator", IOHandler.currentFile());

	IOHandler.writeAttribute("Type",            hdf5String,        &myHead.pType.at(0),  gid);
	IOHandler.writeAttribute("nSites",          H5T_NATIVE_INT,    &myHead.nSize,        gid);
	IOHandler.writeAttribute("Precision",       hdf5String,        &myHead.prec.at(0),   gid);
	IOHandler.writeAttribute("Walltime",        H5T_NATIVE_DOUBLE, &myHead.wallTime,     gid);

	H5Gclose (gid);

	return 0;
}

	Hdf5ReadWriter::Hdf5ReadWriter (iParms myParms) : outputDir(myParms.outputDir), outputName(myParms.outputName), fIndex(myParms.fIndex) {

	char cIndex[8];

	if (fIndex < 0) {
		fIndex = 0;
	}

	sprintf (cIndex, ".%05d\0", fIndex);
	std::string	fileName(outputDir + outputName + std::string(cIndex));

        /*      Create the file and release the plist   */
	if ((file_id = H5Fcreate (fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0)
	{
		printf ("Error opening file %s\n", fileName.c_str());
		return;
	}
}

	Hdf5ReadWriter::~Hdf5ReadWriter () {

	//H5Fclose(file_id);
}

herr_t	Hdf5ReadWriter::nextFile	(int jump) {
	//printf("Im closing file %d\n",file_id);
	//H5Fclose(file_id);

	fIndex += jump + 1;

	char cIndex[8];

	sprintf (cIndex, ".%05d\0", fIndex);
	std::string	fileName(outputDir + outputName + std::string(cIndex));

	if ((file_id = H5Fcreate (fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0)
	{
		printf ("Error opening file %s\n", fileName.c_str());
		return	-1;
	}

	return	0;
}

herr_t  Hdf5ReadWriter::readAttribute	(std::string attName, hid_t h5Type, void *data, hid_t cid)
{
        hid_t   attr;
        herr_t  status;

        if ((attr   = H5Aopen_by_name (cid, ".", attName.c_str(), H5P_DEFAULT, H5P_DEFAULT)) < 0)
                printf ("Error opening attribute %s\n", attName.c_str());
        if ((status = H5Aread (attr, h5Type, data)) < 0)
                printf ("Error reading attribute %s\n", attName.c_str());
        status = H5Aclose(attr);

        return  status;
}

herr_t  Hdf5ReadWriter::readAttribute	(std::string attName, hid_t h5Type, void *data) {
	readAttribute	(attName, h5Type, data, file_id);
	return 0;
}

herr_t  Hdf5ReadWriter::writeAttribute	(std::string attName, hid_t h5Type, void *data, hid_t cid)
{
        hid_t   attr, attr_id;
        herr_t  status;

        attr_id = H5Screate(H5S_SCALAR);
        if ((attr   = H5Acreate2 (cid, attName.c_str(), h5Type, attr_id, H5P_DEFAULT, H5P_DEFAULT)) < 0)
                printf ("Error creating attribute %s\n", attName.c_str());
        if ((status = H5Awrite (attr, h5Type, data)) < 0)
                printf ("Error writing attribute %s to file\n", attName.c_str());
        H5Sclose(attr_id);
        status = H5Aclose(attr);

        return  status;
}

herr_t  Hdf5ReadWriter::writeAttribute	(std::string attName, hid_t h5Type, void *data) {
	writeAttribute	(attName, h5Type, data, file_id);
	return 0;
}

herr_t	Hdf5ReadWriter::writeData (std::string dataName, hid_t h5Type, void *aData, size_t aSize, hid_t cid)
{
	hid_t   base_id, dataSpace, sSpace, dataSet;
	hsize_t dims[1] = { aSize };
	

	base_id = createGroup (dataName, cid);

	/*      Create dataset  */
	dataSpace = H5Screate_simple(1, dims, NULL);
	dataSet   = H5Dcreate(base_id, "data", h5Type, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	sSpace    = H5Dget_space (dataSet);


	/*      Write spectrum data     */
	if (H5Dwrite(dataSet, h5Type, dataSpace, sSpace, H5P_DEFAULT, aData) < 0) {
		printf ("Error writing %zu elements to dataset\n", aSize);
		return	-1;
	}

	/*      Close everything        */
	H5Sclose (sSpace);
	H5Dclose (dataSet);
	H5Sclose (dataSpace);
	H5Gclose (base_id);
	return 0;
}

herr_t	Hdf5ReadWriter::writeData (std::string dataName, hid_t h5Type, void *aData, size_t aSize) {
	writeData (dataName, h5Type, aData, aSize, file_id);
	return 0;
}

herr_t	Hdf5ReadWriter::writeConf (Cosmos *bck, Axiton *field) {

	Hdf5Header	myHeader;

	myHeader.lSize    = bck->CosmosSize();
	myHeader.nQcd     = bck->InitParms().nQcd;
	myHeader.z        = field->z();

	switch (bck->Expansion()) {
		case	Minkowski:

		myHeader.R   = field->R<Minkowski>();
		myHeader.exp = "Minkowski";
		break;

		case	Radiation:

		myHeader.R   = field->R<Radiation>();
		myHeader.exp = "Radiation";
		break;
	}

	switch (bck->InitParms().cType) {
		case	IcFlat:
		myHeader.icType	= "Flat";
		break;

		case	IcSinc2:
		myHeader.icType	= "Sinc2";
		break;

		case	IcOvr2:
		myHeader.icType	= "1/r2";
		break;

		case	IcOvX:
		myHeader.icType	= "1/(1+r)";
		break;

		case	IcGen:
		myHeader.icType	= "Gen";
		break;

		case 	IcNone:
		myHeader.icType = "None";
		break;
	}

	myHeader.icParm1  = bck->InitParms().parm1;
	myHeader.icParm2  = bck->InitParms().parm2;
	myHeader.pType    = "RKN4"; // Use switch when the other propagators are ready
	myHeader.nSize    = bck->CosmosLatt();
	myHeader.massA    = sqrt(bck->AxionMassSq(myHeader.R));

	switch (field->Precision()) 
	{
		case	SinglePrecision:
		myHeader.prec	= "Single";
		break;

		case	DoublePrecision:
		myHeader.prec	= "Double";
		break;
	}

	myHeader.wallTime = bck->InitParms().wTime;


	writeHeader (*this, myHeader);


	auto h5Type = H5T_NATIVE_CHAR;

	switch (field->Precision()) {
		case SinglePrecision:
		h5Type = H5T_NATIVE_FLOAT;
		break;

		case DoublePrecision:
		h5Type = H5T_NATIVE_DOUBLE;
		break;
	}


	if (field->fieldStatus() & FieldGpu == 0)
		field->transferField (FieldBaseDev, HostToDevice);

	writeData ("Field", h5Type, field->fieldCpu(), field->Size());
	writeData ("Dev",   h5Type, field->devCpu(),   field->Size());
	
	//printf("Im closing file %d \n",file_id);
	H5Fclose(file_id);

	return 0;
}

herr_t	Hdf5ReadWriter::readConf  (Cosmos *bck, Axiton *field) { return	0; }
