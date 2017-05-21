#include "PLYReader.h"

PLYReader::PLYReader( char* fname, float scl )
{
	if ( ! ( fin = fopen( fname, "rb" ) ) )
	{
		printf("Unable to open file %s\n", fname);
		return;
	};

	newVertices = new vector<float*>;
	newFaces = new vector<int*>;
	vtxFlag = new vector<bool>;
	vtxExtraDataIndex = new vector<int>;	// extra data index for each non-duplicate vertex
	vtxNewExtraData = new vector<vector<float>>;
	modelname = new char[strlen(fname)];
	strcpy(modelname, fname);

	countDataSize = 1;

	// Parse header
	int vmode = 0, lines = 0 ;
	this->numTrians = this->numVerts = this->extraVertex = this->extraTypeNumber = 0 ;
	char seps[] = " ,\t\n\r ";
	seps[5] = 10 ;
	while ( 1 )
	{
		char header[1024] ;
		// Read in a line (unless the line length is larger than 1024 characters)
		fgets( header, 1024, fin ) ;
		// printf("%s\n", header) ;

		// split the line into
		char* token = strtok( header, seps ) ;
		// token[ strlen(token) - 1 ] = '\0' ;
		// int comp = strcmp ( token, "end_header" ) ;
		// printf("%s %d\n", token, comp ) ;
		// char c = getchar() ;

		////////////////////////////////////////////////////////////
		/* Read headers */
		if ( ! strcmp ( token, "end_header" ) )
		{
			// Header finished
			// printf("End encountered!") ;
			break ;
		} 
		else if ( ! strcmp( token, "format" ) ) 
		{// Format specification				
			token = strtok ( NULL, seps ) ;
			if ( !strcmp ( token, "ascii" ) )
			{// ascii type
				this->mode = 0 ;
			}
			else if ( ! strcmp ( token, "binary_big_endian") )
			{// binary_big_endian type 
				this->mode = 2 ;
			}
			else
			{// binary_little_endian type
				this->mode = 1 ;
			}
		} 
		else if ( ! strcmp( token, "element") )
		{// element specification				
			vmode = 0 ;
			lines = 0 ;
			token = strtok( NULL, seps ) ;
			if ( !strcmp( token, "vertex" ) )
			{// Vertex count					
				token = strtok( NULL, seps ) ;
				sscanf( token, "%d", &(this->numVerts) ) ;

				// Enter vertex mode
				vmode = 1 ;
			} 
			else if ( !strcmp( token, "face" ) )
			{// Face count
				token = strtok( NULL, seps ) ;
				sscanf( token, "%d", &(this->numTrians) ) ;

				// Enter face mode
				vmode = 2 ;
			}
		} 
		else if ( ! strcmp( token, "property" ) ) 
		{
			switch ( vmode )
			{
			case 1: // Vertex mode
				if ( lines >= 3 )
				{   // This means there are more than x,y,z coordinates
					// Extra storage for each vertex
					token = strtok( NULL, seps ) ;
					for ( int i = 0 ; i < totalTypes ; i ++ )
					{
						if ( ! strcmp( token, typeName[i] ) )
						{
							this->extraVertex += typeSize[i] ;
							break ;
						}
					}
					extraTypeNumber++;
					token = strtok( NULL, seps ) ;
					int m = strlen(token);
					char *charTemp = new char[m+1];
					strcpy(charTemp, token);
					// NOTE: You MUST set the last char as '\0', or the free-memory step of this parameter will be wrong.
					charTemp[m] = '\0';
					vtxExtraNames.push_back(charTemp);
				}
				lines ++ ;
				break ;

			case 2: // Face mode
				token = strtok( NULL, seps ) ;
				if (! strcmp( token, "list" ) )
				{
					token = strtok( NULL, seps ) ;
					// Get count data type and size.
					// This is used for reading faces in binary mode
					for ( int i = 0 ; i < totalTypes ; i ++ )
					{
						if ( ! strcmp( token, typeName[i] ) )
						{
							countDataSize = typeSize[i];
							break ;
						}
					}
				}
				break ;
			}
		}
	}
	//printf("Header info: %d VERTICES, %d faces, %d extra bytes per vertex. At %d\n", this->numVerts, this->numTrians, this->extraVertex, ftell( fin ) );
	this->vertArray = new float [ 3 * this->numVerts ] ;
	this->faceArray = new int [ 3 * this->numTrians ] ;
	
	// Assign memory for extra data for each vertex
	if (this->extraTypeNumber != 0)
	{
		vtxExtraData.resize(extraTypeNumber);
		for (int i = 0; i != extraTypeNumber; ++i)
		{
			// Initialize all vertex extra data as 1.0
			vtxExtraData[i].resize(numVerts, 0.0);
		}
	}

	////////////////////////////////////////////////////////////
	/* Read VERTICES and faces*/
	int i;
	char temp[1024] ;
	if ( mode > 0 )
	{
		// Read VERTICES in binary mode
		for ( i = 0 ; i < this->numVerts ; i ++ )
		{
			int haveread = 0 ;
			if ( ( haveread = fread( &(vertArray[ 3 * i ]), sizeof( float ), 3, fin ) ) != 3 )
			{
				if( feof( fin ) )      
				{
					printf( "Read error %d ", ftell( fin ) );
				}

				printf("%d\n", haveread ) ;
				exit(0) ;
			}
			if ( mode == 2 )
			{
				// Flip bits
				for ( int j = 0 ; j < 3 ; j ++ )
				{
					flipBits32( & (vertArray[3*i + j]) ) ;
				}
			}

			//printf("%f,%f,%f\n", vertArray[3*i],vertArray[3*i+1], vertArray[3*i+2]);

			// Read extra data for each vertex
			if (extraTypeNumber != 0)
			{
				for (int j = 0; j != extraTypeNumber; ++j)
				{
					fread( &(vtxExtraData[j][i]), 1, this->extraVertex, fin ) ;
				}
			}
			// Calculate the maximum and minimum coordinates
			if ( i == 0 )
			{
				for ( int j = 0 ; j < 3 ; j ++ )
				{
					rawmin[j] = rawmax[j] = vertArray[3 * i + j] ;
				}
			}
			else
			{
				for ( int k = 0 ; k < 3 ; k ++ )
				{
					if ( rawmin[k] > vertArray[3 * i + k] )
					{
						rawmin[k] = vertArray[3 * i + k] ;
					}
					if ( rawmax[k] < vertArray[3 * i + k] )
					{
						rawmax[k] = vertArray[3 * i + k] ;
					}
				}
			}

		}

		// Read faces in binary mode
		for ( i = 0 ; i < this->numTrians ; i ++ )
		{
			int haveread = 0;
			int numVtxPerFace;
			if ( ( haveread = fread(&numVtxPerFace, countDataSize, 1, fin )) != 1)
			{// Read count data for a face
				if( feof( fin ) )      
				{
					printf( "Read error %d ", ftell( fin ) );
				}
				printf("%d\n", haveread ) ;
				exit(0) ;
			}

			if ( ( haveread = fread(&(faceArray[3*i]), sizeof(int), 3, fin )) != 3 )
			{// Read face data for a face
				if( feof( fin ) )      
				{
					printf( "Read error %d ", ftell( fin ) );
				}
				printf("%d\n", haveread ) ;
				exit(0) ;
			}

			//printf("%d,%d,%d\n", faceArray[3*i],faceArray[3*i+1], faceArray[3*i+2]);
		}

	}
	else
	{
		// Read VERTICES in Ascii mode
		for ( i = 0 ; i < this->numVerts ; i ++ )
		{
			// Read coordinates x,y,z
			fgets ( temp, 1024, fin ) ;
			char seps[] = " ,\t\n";
			char* token = strtok( temp, seps ) ;

			sscanf( token, "%f", &(vertArray[ 3 * i ] ) ) ;
			token = strtok( NULL, seps ) ;
			sscanf( token, "%f", &(vertArray[ 3 * i + 1 ] ) ) ;
			token = strtok( NULL, seps ) ;
			sscanf( token, "%f", &(vertArray[ 3 * i + 2 ] ) ) ;

			// Read extra data for each vertex
			if (extraTypeNumber != 0)
			{
				for (int j = 0; j != extraTypeNumber; ++j)
				{
					token = strtok( NULL, seps ) ;
					sscanf( token, "%f", &(vtxExtraData[j][i]) ) ;
				}
			}

			// Compute the maximum and minimum coordinates for each dimension
			if ( i == 0 )
			{
				for ( int j = 0 ; j < 3 ; j ++ )
				{
					rawmin[j] = rawmax[j] = vertArray[3 * i + j] ;
				}
			}
			else
			{
				for ( int k = 0 ; k < 3 ; k ++ )
				{
					if ( rawmin[k] > vertArray[3 * i + k] )
					{
						rawmin[k] = vertArray[3 * i + k] ;
					}
					if ( rawmax[k] < vertArray[3 * i + k] )
					{
						rawmax[k] = vertArray[3 * i + k] ;
					}
				}
			}
		}

		// Read faces in ASCII mode
		for (i = 0; i < this->numTrians; ++i)
		{
			fgets ( temp, 1024, fin ) ;
			char seps[] = " ,\t\n";
			char* token = strtok( temp, seps ) ;

			int numVtxPerFace = 0;
			sscanf( token, "%d", &numVtxPerFace) ;
			if (numVtxPerFace > 3)
			{
				printf("This model contains polygons with more than 3 edges. Quitting\n");
				break;
			}

			token = strtok( NULL, seps ) ;
			sscanf( token, "%d", &(faceArray[ 3 * i ] ) ) ;
			token = strtok( NULL, seps ) ;
			sscanf( token, "%d", &(faceArray[ 3 * i + 1 ] ) ) ;
			token = strtok( NULL, seps ) ;
			sscanf( token, "%d", &(faceArray[ 3 * i + 2 ] ) ) ;
		}
	}

	// determine the maximum and minimum coordinates
	maxsize = rawmax[0] - rawmin[0] ;
	if ( rawmax[1] - rawmin[1] > maxsize )
	{
		maxsize = rawmax[1] - rawmin[1] ;
	}
	if ( rawmax[2] - rawmin[2] > maxsize )
	{
		maxsize = rawmax[2] - rawmin[2] ;
	}

	for ( i = 0 ; i < 3 ; i ++ )
	{
		min[i] = ( rawmax[i] + rawmin[i] ) / 2 - maxsize / 2 ;
		max[i] = ( rawmax[i] + rawmin[i] ) / 2 + maxsize / 2 ;
	}

	// printf( "Low corner: %f %f %f   High corner %f %f %f \n", min[0], min[1], min[2], max[0], max[1], max[2] );
	// printf( "Maxdif: %f \n", maxsize );

	// Scale
	for ( i = 0 ; i < 3 ; i ++ )
	{
		min[i] -= maxsize * ( 1 / scl - 1 ) / 2 ;
	}
	maxsize *= ( 1 / scl ) ;

	// Reset 
	this->offset = ftell( fin ) ;
	curtrian = 0 ;

	computeCenter();
}

void PLYReader::computeCenter()
{
	ctr[0] = ctr[1] = ctr[2] = 0;
	for (int i = 0; i != numVerts; ++i)
	{
		float *v = getVertex(i);
		for (int j = 0; j != 3; ++j)
		{
			ctr[j] += v[j];
		}
	}
	for (int i = 0; i != 3; ++i)
	{
		ctr[i] /= numVerts;
	}
}

void PLYReader::DeleteDuplicateVertices()
{
	//////////////////////////////////////////////////////////////////////////
	// Create a K-D tree and add all VERTICES into the tree
	struct kdtree *kd = kd_create(3);
	struct kdres *presults;
	float *v;
	double vtx[3];
	int *data;
	for (int i = 0; i < numVerts; ++i)
	{
		v = getVertex(i);
		vtx[0] = (double)v[0];
		vtx[1] = (double)v[1];
		vtx[2] = (double)v[2];
		data = new int;
		*data = i;
		kd_insert3(kd, vtx[0],vtx[1],vtx[2],data) ;
	}
	//////////////////////////////////////////////////////////////////////////
	// Get the nearest VERTICES of each vertex
	const double radius = 0.0000001;
	double pos[3];
	vector<int> *vecVtxIndex = new vector<int>;	// nearest vertex indexes
	for (int i = 0; i < numVerts; ++i)
	{// initialize as -1
		vecVtxIndex->push_back(-1);
	}
	for (int i = 0; i < numVerts; ++i)
	{
		if (vecVtxIndex->at(i) == -1)
		{// the nearest vertex still needs to be determined
			v = getVertex(i);
			vtx[0] = (double)v[0];
			vtx[1] = (double)v[1];
			vtx[2] = (double)v[2];
			presults = kd_nearest_range( kd, vtx, radius);

			while( !kd_res_end( presults ) )
			{// Process all nearest VERTICES
				int *ptr_nearestVtx = (int*)kd_res_item( presults, pos );
				int nearestVtx = *ptr_nearestVtx;
				vecVtxIndex->at(nearestVtx) = i;

				/* go to the next entry */
				kd_res_next( presults );
			}
			kd_res_free( presults );
		}
	}

	//delete data;
	kd_free( kd );

	//////////////////////////////////////////////////////////////////////////
	/* Get non-duplicate VERTICES */
	
	// The ith value in this vector "vecSmallNonDupCount" is the number of non-duplicate VERTICES with index smaller 
	// than i for the ith vertex. This kind of value is used to update the new vertex indexes in
	// the new face connectivity. For instance, for the 8th vertex, if 5 VERTICES are non-duplicate VERTICES with
	// indexes smaller than 8 (such as 0, 1, 2, 4, 5), then the 8th value in "vecSmallNonDupCount" is 5. This
	// means that the new index for the 8th vertex in the new face connectivity is 5 instead of 8.
	vector<int> *vecSmallNonDupCount = new vector<int>;	
	int count = 0;
	for (int i = 0; i < numVerts; ++i)
	{
		if (vecVtxIndex->at(i) == i)
		{// This index is non-duplicate
			v = getVertex(i);
			float *v_temp = new float[3];
			v_temp[0] = v[0];
			v_temp[1] = v[1];
			v_temp[2] = v[2];
			newVertices->push_back(v_temp);
			vtxExtraDataIndex->push_back(i);
			vecSmallNonDupCount->push_back(count);	// save updated index for face
			count++;	// real index plus 1
		}
		else
		{// This vertex is duplicate, so set its index as -1
			vecSmallNonDupCount->push_back(-1);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Update face indexes and delete degenerate faces
	//newFaces = new vector<int*>;
	for (int i = 0; i < numTrians; ++i)
	{
		int *t = getTriangle(i);
		int *t_temp = new int[3];
		for (int j = 0; j < 3; ++j)
		{
			if (vecVtxIndex->at(t[j]) != t[j])
			{// If a vertex is duplicate, updated its index in the face
				t_temp[j] = vecVtxIndex->at(t[j]);
			}
			else
			{
				t_temp[j] = t[j];
			}
			// Get Newly updated vertex index in this face
			t_temp[j] = vecSmallNonDupCount->at(t_temp[j]);
		}
		// Check the degenerate faces
		if ((t_temp[0] != t_temp[1]) && (t_temp[0] != t_temp[2]) && (t_temp[1] != t_temp[2]))
		{// if the face is not a degenerate face, save it
			newFaces->push_back(t_temp);
		}		
	}
	vecVtxIndex->clear();
	vecSmallNonDupCount->clear();
}

void PLYReader::writeNewPLY(char *fname)
{
	// Note: "wb" denotes to write a binary file, 
	//	"w" denotes to write a ascii file
	FILE *fout;
	
	if (mode == 0)
	{
		fout = fopen( fname, "w" );
	}
	else
	{
		fout = fopen( fname, "wb" );
	}
	if ( fout == NULL )
	{
		printf("Unable to create file %s\n", fname);
		return;
	}
	////////////////////////////////////////////////////////////
	// Write headers
	// Ply
	fprintf( fout, "ply\n" ) ;
	if (mode == 1)
	{
		fprintf( fout, "format binary_little_endian 1.0\n" ) ;
	}
	else if (mode == 2)
	{
		fprintf( fout, "format binary_big_endian 1.0\n" ) ;
	}
	else
	{
		fprintf( fout, "format ascii 1.0\n" ) ;
	}

	// vertex properties
	fprintf( fout, "element vertex %d\n", (int)newVertices->size() ) ;
	fprintf( fout, "property float x\n" ) ;
	fprintf( fout, "property float y\n" ) ;
	fprintf( fout, "property float z\n" ) ;

	if (this->extraTypeNumber != 0)
	{
		for (int i = 0; i != this->extraTypeNumber; ++i)
		{
			fprintf( fout, "property float %s\n", vtxExtraNames[i]) ;
		}
	}

	// face properties
	fprintf( fout, "element face %d\n", (int)newFaces->size() ) ;
	fprintf( fout, "property list uchar int vertex_indices\n" ) ;

	// End
	fprintf( fout, "end_header\n" ) ;

	////////////////////////////////////////////////////////////
	// Write VERTICES
	float *v;
	for (int i = 0; i != newVertices->size(); ++i)
	{
		v = newVertices->at(i);
		if (mode == 2)
		{
			flipBits32( &(v[i]) ) ;
		}
		if (mode == 0)
		{
			if (this->extraTypeNumber == 0)
			{
				fprintf( fout, "%f %f %f\n",v[0], v[1], v[2]) ;
			}
			else
			{
				fprintf( fout, "%f %f %f",v[0], v[1], v[2]) ;
				for (int j = 0; j != this->extraTypeNumber; ++j)
				{
					fprintf( fout, " %f", vtxExtraData[j][vtxExtraDataIndex->at(i)]);
				}
				fprintf( fout, "\n");
			}
		}
		else
		{
			if (this->extraTypeNumber == 0)
			{
				fwrite( v, sizeof(float), 3, fout ) ;
			}
			else
			{
				fwrite( v, sizeof(float), 3, fout ) ;
				for (int j = 0; j != this->extraTypeNumber; ++j)
				{
					fwrite( &(vtxExtraData[j][vtxExtraDataIndex->at(i)]), this->extraVertex, 1, fout ) ;
				}
			}
			
		}			
	}
	////////////////////////////////////////////////////////////
	// Write Faces
	int f[3];
	int *t;
	int cnum = 3;

	for ( int i = 0 ; i != newFaces->size(); i++ )
	{
		t = newFaces->at(i);
		for (int j = 0; j < 3; ++j)
		{
			f[j] = t[j];
		}
		if (mode == 0)
		{
			fprintf( fout, "%u ",cnum);
			fprintf( fout, "%d %d %d\n",f[0], f[1], f[2]);
		}
		else
		{
			PLYWriter::writeFace(fout, cnum, f);
		}
	}
	fclose(fout);
}


void PLYReader::writePLYWithExtraData(char *fname, vector<char*> vecExtraNames, vector<vector<float>> vecExtraData)
{
	// Note: "wb" denotes to write a binary file, 
	//	"w" denotes to write a ascii file
	FILE *fout;
	if (mode == 0)
	{
		fout = fopen( fname, "w" );
	}
	else
	{
		fout = fopen( fname, "wb" );
	}	
	
	if ( fout == NULL )
	{
		printf("Unable to create file %s\n", fname);
		return;
	}
	////////////////////////////////////////////////////////////
	// Write headers
	// Ply
	fprintf( fout, "ply\n" ) ;
	if (mode == 1)
	{
		fprintf( fout, "format binary_little_endian 1.0\n" ) ;
	}
	else if (mode == 2)
	{
		fprintf( fout, "format binary_big_endian 1.0\n" ) ;
	}
	else
	{
		fprintf( fout, "format ascii 1.0\n" ) ;
	}

	// vertex properties
	fprintf( fout, "element vertex %d\n", numVerts ) ;
	fprintf( fout, "property float x\n" ) ;
	fprintf( fout, "property float y\n" ) ;
	fprintf( fout, "property float z\n" ) ;
	for (int i = 0; i != vecExtraNames.size(); ++i)
	{
		fprintf( fout, "property float %s\n", vecExtraNames[i]) ;
	}

	// face properties
	fprintf( fout, "element face %d\n", numTrians ) ;
	fprintf( fout, "property list uchar int vertex_indices\n" ) ;

	// End
	fprintf( fout, "end_header\n" ) ;

	////////////////////////////////////////////////////////////
	// Write VERTICES
	float *v;
	for (int i = 0; i != numVerts; ++i)
	{
		v = getVertex(i);
		if (mode == 2)
		{
			flipBits32( &(v[i]) ) ;
		}
		if (mode == 0)
		{
			fprintf( fout, "%f %f %f",v[0], v[1], v[2]) ;
			for (int j = 0; j != vecExtraData.size(); ++j)
			{
				fprintf( fout, " %f", vecExtraData[j][i]);
			}
			fprintf( fout, "\n");
		}
		else
		{
			fwrite( v, sizeof(float), 3, fout ) ;
			for (int j = 0; j != vecExtraData.size(); ++j)
			{
				fwrite( &(vecExtraData[j][i]), 4, 1, fout ) ;
			}
		}			
	}
	////////////////////////////////////////////////////////////
	// Write Faces
	int f[3];
	int *t;
	int cnum = 3;

	for ( int i = 0 ; i != numTrians; i++ )
	{
		t = getTriangle(i);
		for (int j = 0; j < 3; ++j)
		{
			f[j] = t[j];
		}
		if (mode == 0)
		{
			fprintf( fout, "%u ",cnum);
			fprintf( fout, "%d %d %d\n",f[0], f[1], f[2]);
		}
		else
		{
			unsigned char cnumval = cnum ;
			fwrite( &cnumval, sizeof( unsigned char ), 1, fout ) ;
			fwrite( f, sizeof(int), cnum, fout ) ;
		}
	}
	fclose(fout);
}


void PLYReader::writeNewPLYWithExtraData(char *fname)
{
	// Note: "wb" denotes to write a binary file, 
	//	"w" denotes to write a ascii file
	FILE *fout;
	//mode = 0;
	if (mode == 0)
	{
		fout = fopen( fname, "w" );
	}
	else
	{
		fout = fopen( fname, "wb" );
	}	

	if ( fout == NULL )
	{
		printf("Unable to create file %s\n", fname);
		return;
	}
	////////////////////////////////////////////////////////////
	// Write headers
	// Ply
	fprintf( fout, "ply\n" ) ;
	if (mode == 1)
	{
		fprintf( fout, "format binary_little_endian 1.0\n" ) ;
	}
	else if (mode == 2)
	{
		fprintf( fout, "format binary_big_endian 1.0\n" ) ;
	}
	else
	{
		fprintf( fout, "format ascii 1.0\n" ) ;
	}

	// vertex properties
	fprintf( fout, "element vertex %d\n", newVertices->size() ) ;
	fprintf( fout, "property float x\n" ) ;
	fprintf( fout, "property float y\n" ) ;
	fprintf( fout, "property float z\n" ) ;
	for (int i = 0; i != extraTypeNumber; ++i)
	{
		fprintf( fout, "property float %s\n", vtxExtraNames[i]) ;
	}
	// face properties
	fprintf( fout, "element face %d\n", newFaces->size() ) ;
	fprintf( fout, "property list uchar int vertex_indices\n" ) ;

	// End
	fprintf( fout, "end_header\n" ) ;

	////////////////////////////////////////////////////////////
	// Write VERTICES
	float *v;
	for (int i = 0; i != newVertices->size(); ++i)
	{
		v = newVertices->at(i);
		if (mode == 2)
		{
			flipBits32( &(v[i]) ) ;
		}
		if (mode == 0)
		{
			fprintf( fout, "%f %f %f",v[0], v[1], v[2]) ;
			for (int j = 0; j != extraTypeNumber; ++j)
			{
				fprintf( fout, " %f", (vtxNewExtraData->at(j))[i]);
			}
			fprintf( fout, "\n");
		}
		else
		{
			fwrite( v, sizeof(float), 3, fout ) ;
			for (int j = 0; j != vtxNewExtraData->size(); ++j)
			{
				fwrite( &((vtxNewExtraData->at(j))[i]), 4, 1, fout ) ;
			}
		}			
	}
	////////////////////////////////////////////////////////////
	// Write Faces
	int f[3];
	int *t;
	int cnum = 3;

	for ( int i = 0 ; i != newFaces->size(); i++ )
	{
		t = newFaces->at(i);
		for (int j = 0; j < 3; ++j)
		{
			f[j] = t[j];
		}
		if (mode == 0)
		{
			fprintf( fout, "%u ",cnum);
			fprintf( fout, "%d %d %d\n",f[0], f[1], f[2]);
		}
		else
		{
			unsigned char cnumval = cnum ;
			fwrite( &cnumval, sizeof( unsigned char ), 1, fout ) ;
			fwrite( f, sizeof(int), cnum, fout ) ;
		}
	}
	fclose(fout);
}


void PLYReader::writeOriginalPLY(char *fname)
{
	FILE *fout = NULL;
	if (mode == 0)
	{// "w" denotes to write a ascii file
		fout = fopen( fname, "w" );
	}
	else
	{// "wb" denotes to write a binary file
		fout = fopen( fname, "wb" );
	}
	//FILE *fout = fopen( fname, "wb" );
	if ( fout == NULL )
	{
		printf("Unable to create file %s\n", fname);
		return;
	}
	////////////////////////////////////////////////////////////
	// Write headers
	fprintf( fout, "ply\n" ) ;
	if (mode == 1)
	{
		fprintf( fout, "format binary_little_endian 1.0\n" ) ;
	}
	else if (mode == 2)
	{
		fprintf( fout, "format binary_big_endian 1.0\n" ) ;
	}
	else
	{
		fprintf( fout, "format ascii 1.0\n" ) ;
	}
	// vertex properties
	fprintf( fout, "element vertex %d\n", numVerts) ;
	fprintf( fout, "property float x\n" ) ;
	fprintf( fout, "property float y\n" ) ;
	fprintf( fout, "property float z\n" ) ;
	// face properties
	fprintf( fout, "element face %d\n", numTrians) ;
	fprintf( fout, "property list uchar int vertex_indices\n" ) ;
	// End
	fprintf( fout, "end_header\n" ) ;

	////////////////////////////////////////////////////////////
	// Write VERTICES
	float *v;
	for (int i = 0; i != numVerts; ++i)
	{
		v = getVertex(i);
		if (mode == 2)
		{
			flipBits32( &(v[i]) ) ;
		}
		if (mode == 0)
		{
			fprintf( fout, "%f %f %f\n",v[0], v[1], v[2]) ;
		}
		else
		{
			fwrite( v, sizeof(float), 3, fout ) ;
		}			
	}
	////////////////////////////////////////////////////////////
	// Write Faces
	int f[3];
	int *t;
	int cnum = 3;

	for ( int i = 0 ; i != numTrians; i++ )
	{
		t = getTriangle(i);
		for (int j = 0; j < 3; ++j)
		{
			f[j] = t[j];
		}
		if (mode == 0)
		{
			fprintf( fout, "%u ",cnum);
			fprintf( fout, "%d %d %d\n",f[0], f[1], f[2]);
		}
		else
		{
			unsigned char cnumval = cnum ;
			fwrite( &cnumval, sizeof( unsigned char ), 1, fout ) ;
			fwrite( f, sizeof(int), cnum, fout ) ;
		}
	}
	fclose(fout);
}

void PLYReader::flipBits32( void *x )
{
	unsigned char *temp = (unsigned char *)x;
	unsigned char swap;

	swap = temp [ 0 ];
	temp [ 0 ] = temp [ 3 ];
	temp [ 3 ] = swap;

	swap = temp [ 1 ];
	temp [ 1 ] = temp [ 2 ];
	temp [ 2 ] = swap;
}

void PLYReader::clearNewVertexFaceVector()
{
	if (!newVertices->empty())
	{
		for (int i = 0; i != newVertices->size(); ++i)
		{
			delete (*newVertices)[i];
		}
		newVertices->clear();
		ClearVector(*newVertices);
	}
	if (!newFaces->empty())
	{
		for (int i = 0; i != newFaces->size(); ++i)
		{
			delete (*newFaces)[i];
		}
		newFaces->clear();
		ClearVector(*newFaces);
	}
	if (!vtxNewExtraData->empty())
	{
		for (int i = 0; i != vtxNewExtraData->size(); ++i)
		{
			vtxNewExtraData->at(i).clear();
			ClearVector(vtxNewExtraData->at(i));
		}
		vtxNewExtraData->clear();
		ClearVector(*vtxNewExtraData);
	}
}

PLYReader::~PLYReader()
{
	delete []vertArray;
	delete []faceArray;
	clearNewVertexFaceVector();
	if (!vtxFlag->empty())
	{
		vtxFlag->clear();
		ClearVector(*vtxFlag);
	}
	for (int i = 0; i != vtxExtraData.size();++i)
	{
		vtxExtraData[i].clear();
		ClearVector(vtxExtraData[i]);
	}
	vtxExtraData.clear();
	ClearVector(vtxExtraData);

	for (int i = 0; i != vtxExtraNames.size(); ++i)
	{
		delete[] vtxExtraNames[i];
		vtxExtraNames[i] = NULL;
	}
	vtxExtraNames.clear();
	ClearVector(vtxExtraNames);
	vtxExtraDataIndex->clear();
	ClearVector(*vtxExtraDataIndex);
	if (modelname != NULL)
	{
		delete modelname;
	}
}

void PLYReader::RemoveLargestDifferenceFrom(PLYReader *inputPly)
{
	vector<int> *vecFaceIndex = new vector<int>;

	cout<<endl<<"Getting difference between two meshes..."<<endl;
	GetDifferenceWithAnotherMesh(inputPly,vecFaceIndex);
	cout<<"  Total different face number:"<<vecFaceIndex->size()<<endl;

	vector<int> *vecResultVtxIndex = new vector<int>;
	
	const int num = 1;
	cout<<"Getting "<<num<<" largest connected components of the difference data..."<<endl;
	GetConnectedComponents(vecFaceIndex, 0, num, vecResultVtxIndex);
	
	cout<<"Deleting "<<num<<" largest connected component from the difference data..."<<endl;
	DeleteVertices(vecResultVtxIndex);
	vecFaceIndex->clear();
	vecResultVtxIndex->clear();
	
}

void PLYReader::GetDifferenceWithAnotherMesh(PLYReader *inputPly, vector<int> *vecResultFaceIndex)
{
	// Create a K-D tree based on all VERTICES from the original model
	const double radius = 0.08;
	struct kdtree *kd = kd_create(3);
	struct kdres *presults;
	float *v;
	double vtx[3];
	int *data;
	for (int i = 0; i < numVerts; ++i)
	{
		v = getVertex(i);
		vtx[0] = (double)v[0];
		vtx[1] = (double)v[1];
		vtx[2] = (double)v[2];
		data = new int;
		*data = i;
		kd_insert3(kd, vtx[0],vtx[1],vtx[2],data) ;
	}

	vtxFlag->resize(numVerts, true);
	// Get the nearest VERTICES of each vertex
	double pos[3];
	for (int i = 0; i < inputPly->numVerts; ++i)
	{
		v = inputPly->getVertex(i);
		vtx[0] = (double)v[0];
		vtx[1] = (double)v[1];
		vtx[2] = (double)v[2];
		presults = kd_nearest_range( kd, vtx, radius);

		while( !kd_res_end( presults ))
		{// Process all nearest VERTICES inside the sphere defined by the radius with input vertex as center
			int *ptr_nearestVtx = (int*)kd_res_item( presults, pos );
			int nearestVtx = *ptr_nearestVtx;
			// Put the nearest vertex as false
			if (vtxFlag->at(nearestVtx))
			{
				vtxFlag->at(nearestVtx) = false;
			}			

			/* go to the next entry */
			kd_res_next( presults );
		}
		kd_res_free( presults );
	}

	//delete data;
	kd_free( kd );

	//// Get the different vertex between two meshes
	//vector<int> *vecVtxIndex = new vector<int>;
	//for (int i = 0; i != numVerts; ++i)
	//{
	//	if (vtxFlag->at(i))
	//	{
	//		vecVtxIndex->push_back(i);
	//	}
	//}

	////
	//cout<<"  Number of patching VERTICES:"<<vecVtxIndex->size()<<endl;

	// Get the different faces between two meshes
	for (int i = 0; i != numTrians; ++i)
	{
		int *f = getTriangle(i);
		int flag = true;
		for (int j = 0; j != 3; ++j)
		{
			if (!vtxFlag->at(f[j]))
			{// the triangle contains a vertex that is not what we need
				flag = false;
				break;
			}
		}
		if (flag)
		{// If the triangle contains three VERTICES we need, save it
			vecResultFaceIndex->push_back(i);
		}
	}
}

void PLYReader::computeNearestEuclideanDistance(PLYReader *inputPly, vector<double> *vecDis)
{
	// Create a K-D tree based on all VERTICES from the input model
	const double radius = 0.05;
	struct kdtree *kd = kd_create(3);
	struct kdres *presults;
	float *v;
	double vtx[3];
	double dis = 0.0;
	int *data;
	for (int i = 0; i < inputPly->numVerts; ++i)
	{
		v = inputPly->getVertex(i);
		vtx[0] = (double)v[0];
		vtx[1] = (double)v[1];
		vtx[2] = (double)v[2];
		data = new int;
		*data = i;
		kd_insert3(kd, vtx[0],vtx[1],vtx[2],data) ;
	}

	// Get the nearest VERTICES of each vertex from the original model
	double pos[3];
	for (int i = 0; i < numVerts; ++i)
	{
		v = getVertex(i);
		vtx[0] = (double)v[0];
		vtx[1] = (double)v[1];
		vtx[2] = (double)v[2];
		// Find the nearest vertex of vtx
		presults = kd_nearest (kd, vtx);
		int *ptr_nearestVtx = (int*)kd_res_item( presults, pos );
		dis = sqrt((vtx[0]-pos[0])*(vtx[0]-pos[0])
			+ (vtx[1]-pos[1])*(vtx[1]-pos[1])
			+ (vtx[2]-pos[2])*(vtx[2]-pos[2]));
		vecDis->push_back(dis);
		kd_res_free( presults );
	}

	//delete kd-tree;
	kd_free( kd );
}

void PLYReader::getExtraDataForNewMesh(PLYReader *inputPly, vector<double> *vecDis)
{
	// Compute the maximum, minimum value, mean value and standard deviation for the distances
	double meanVal = 0, stdVal = 0;
	for (int i = 0; i != vecDis->size(); ++i)
	{
		meanVal += vecDis->at(i);
	}
	meanVal = meanVal/vecDis->size();
	for (int i = 0; i != vecDis->size(); ++i)
	{
		stdVal += (vecDis->at(i) - meanVal) * (vecDis->at(i) - meanVal);
	}
	stdVal = stdVal/(vecDis->size()-1);
	stdVal = sqrt(stdVal);

	// Compute a threshold for determining the extra data for each vertex in the input mesh
	const int k_val = 6;
	const double thresholdVal = meanVal + k_val * stdVal;

	// Create a K-D tree based on all VERTICES from the original mesh
	struct kdtree *kd = kd_create(3);
	struct kdres *presults;
	float *v;
	double vtx[3];
	double dis = 0.0;
	int *data;
	for (int i = 0; i < numVerts; ++i)
	{
		v = getVertex(i);
		vtx[0] = (double)v[0];
		vtx[1] = (double)v[1];
		vtx[2] = (double)v[2];
		data = new int;
		*data = i;
		kd_insert3(kd, vtx[0],vtx[1],vtx[2],data) ;
	}

	// Copy some basic extra information from original model to input model
	inputPly->extraVertex = this->extraVertex;
	inputPly->extraTypeNumber = this->extraTypeNumber;
	inputPly->vtxExtraData.resize(this->extraTypeNumber);
	inputPly->vtxExtraNames.resize(this->extraTypeNumber);
	//vtxExtraDataIndex = new vector<int>;

	// Get the nearest VERTICES of each vertex from the new input mesh
	double pos[3];
	for (int i = 0; i != this->extraTypeNumber; ++i)
	{
		inputPly->vtxExtraData[i].resize(inputPly->numVerts, 0);
	}

	for (int i = 0; i < inputPly->numVerts; ++i)
	{
		v = inputPly->getVertex(i);
		vtx[0] = (double)v[0];
		vtx[1] = (double)v[1];
		vtx[2] = (double)v[2];
		// Find the nearest vertex of vtx
		presults = kd_nearest (kd, vtx);
		int *ptr_nearestVtx = (int*)kd_res_item( presults, pos );
		dis = sqrt((vtx[0]-pos[0])*(vtx[0]-pos[0])
			+ (vtx[1]-pos[1])*(vtx[1]-pos[1])
			+ (vtx[2]-pos[2])*(vtx[2]-pos[2]));
		if (dis <= thresholdVal)
		{
			for (int j = 0; j != this->extraTypeNumber; ++j)
			{
				inputPly->vtxExtraData[j][i] = this->vtxExtraData[j][*ptr_nearestVtx];
			}
		}
		else
		{
			for (int j = 0; j != this->extraTypeNumber; ++j)
			{
				inputPly->vtxExtraData[j][i] = 0;
			}
		}
		kd_res_free( presults );
	}

	//delete kd-tree;
	kd_free( kd );
}

void PLYReader::GetConnectedComponents
	(vector<int> *vecFaceIndex, int startIndex, int endIndex, vector<int> *vecResultVtxIndex)
{
	ALGraph *alGraph = new ALGraph(numVerts);
	alGraph->createALGraph(vecFaceIndex, faceArray);

	vector<int> *vecStartNodeID = new vector<int>;
	vector<int> *vecContCompNum = new vector<int>;
	alGraph->GetAllConnectedComponentNumber(vecStartNodeID, vecContCompNum);

	vector<pair<int,int>> *contComp = new vector<pair<int,int>>;
	for (int i = 0; i != vecContCompNum->size(); ++i)
	{
		cout<<"  Connected component start at "<<vecStartNodeID->at(i)
			<<": total number "<<vecContCompNum->at(i)<<endl;
		contComp->push_back(make_pair(vecStartNodeID->at(i),vecContCompNum->at(i)));
	}

	// Sort the connected component number from largest to smallest
	sort(contComp->begin(), contComp->end(), compPair);

	for (int i = startIndex; i != endIndex; ++i)
	{
		cout<<"No. "<<i+1<<" largest connected Component: start at "<<contComp->at(i).first
			<<": total number "<<contComp->at(i).second<<endl;

		alGraph->GetConnectedComponent(contComp->at(i).first);
		for (int j = 0; j != alGraph->vecVtxTraverse->size(); ++j)
		{
			vecResultVtxIndex->push_back(alGraph->vecVtxTraverse->at(j));
		}
	}
}

void PLYReader::GetSmallConnectedComponents(vector<int> *vecFaceIndex, int vertexThreshold, vector<int> *vecResultVtxIndex)
{
	ALGraph *alGraph = new ALGraph(numVerts);
	alGraph->createALGraph(vecFaceIndex, faceArray);

	vector<int> *vecStartNodeID = new vector<int>;
	vector<int> *vecContCompNum = new vector<int>;
	alGraph->GetAllConnectedComponentNumber(vecStartNodeID, vecContCompNum);

	int contCompNum = vecContCompNum->size();
	for (int i = 0; i != contCompNum; ++i)
	{
		if (vecContCompNum->at(i) < vertexThreshold)
		{
			cout<<"No. "<<i+1<<" largest connected Component: start at "<<vecStartNodeID->at(i)
				<<": total number "<<vecContCompNum->at(i)<<endl;
			alGraph->GetConnectedComponent(vecStartNodeID->at(i));
			for (int j = 0; j != alGraph->vecVtxTraverse->size(); ++j)
			{
				vecResultVtxIndex->push_back(alGraph->vecVtxTraverse->at(j));
			}
		}
	}
}

void PLYReader::GetfloatingVertices(vector<int> *vecVtx)
{
	vtxFlag->clear();
	vtxFlag = new vector<bool>;
	vtxFlag->resize(numVerts, false);
	for (int i = 0; i != numTrians; ++i)
	{
		int *f = getTriangle(i);
		for (int j = 0; j != 3; ++j)
		{
			vtxFlag->at(f[j]) = true;
		}
	}
	for (int i = 0; i != numVerts; ++i)
	{
		if (!vtxFlag->at(i))
		{
			vecVtx->push_back(i);
		}
	}
}

void PLYReader::DeleteVertices(vector<int> *vecDeletedVtx)
{
	vtxFlag->clear();
	vtxFlag = new vector<bool>;
	vtxFlag->resize(numVerts, false);
	newVertices->clear();
	newVertices = new vector<float*>;
	newFaces->clear();
	newFaces = new vector<int*>;
	//vtxExtraDataIndex = new vector<int>;

	// Get all VERTICES which needs to be saved
	for (int i = 0; i != vecDeletedVtx->size(); ++i)
	{
		vtxFlag->at(vecDeletedVtx->at(i)) = true;
	}
	int count = 0;
	vector<int> *vecVtxCount = new vector<int>;
	vecVtxCount->resize(numVerts, 0);
	for (int i = 0; i != numVerts; ++i)
	{
		if (!vtxFlag->at(i))
		{
			float *tempVtx = new float[3];
			float *v = getVertex(i);
			for (int j = 0; j != 3; ++j)
			{
				tempVtx[j] = v[j];
			}
			newVertices->push_back(tempVtx);
			vtxExtraDataIndex->push_back(i);
		}
		else
		{
			count++;
		}
		vecVtxCount->at(i) = count;
	}
	
	// Get all triangles which needs to be saved
	for (int i = 0; i != numTrians; ++i)
	{
		int *f = getTriangle(i);
		int flag = true;
		for (int j = 0; j != 3; ++j)
		{
			if (vtxFlag->at(f[j]))
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			int *tempFace = new int[3];
			for (int j = 0; j != 3; ++j)
			{
				tempFace[j] = f[j] - vecVtxCount->at(f[j]);
			}
			newFaces->push_back(tempFace);
		}
	}
	vecVtxCount->clear();
}