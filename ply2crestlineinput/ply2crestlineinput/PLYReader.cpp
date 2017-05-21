#include "PLYReader.h"

PLYReader::~PLYReader()
{
	delete []vertarray;
	delete []faceArray;
}


PLYReader::PLYReader( char* fname )
{
	if ( ! ( fin = fopen( fname, "rb" ) ) )
	{
		printf("Unable to open file %s\n", fname);
		return;
	};

	countDataSize = 1;

	// Parse header
	int vmode = 0, lines = 0 ;
	this->numTrians = this->numVerts = this->extraVertex = 0 ;
	char seps[] = " ,\t\n\r ";
	seps[5] = 10 ;

	////////////////////////////////////////////////////////////
	/* Read headers */
	while ( 1 )
	{
		char header[1024] ;
		// Read in a line (unless the line length is larger than 1024 characters)
		fgets( header, 1024, fin ) ;

		// split the line into
		char* token = strtok( header, seps ) ;

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
	//printf("Header info: %d vertices, %d faces, %d extra bytes per vertex. At %d\n", this->numVerts, this->numTrians, this->extraVertex, ftell( fin ) );
	this->vertarray = new float [ 3 * this->numVerts ] ;
	this->faceArray = new int [ 3 * this->numTrians ] ;

	////////////////////////////////////////////////////////////
	/* Read vertices and faces*/
	int i;
	char temp[1024] ;
	if ( mode > 0 )
	{
		// Read vertices in binary mode
		for ( i = 0 ; i < this->numVerts ; i ++ )
		{
			int haveread = 0 ;
			if ( ( haveread = fread( &(vertarray[ 3 * i ]), sizeof( float ), 3, fin ) ) != 3 )
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
					flipBits32( & (vertarray[3*i + j]) ) ;
				}
			}

			//printf("%f,%f,%f\n", vertarray[3*i],vertarray[3*i+1], vertarray[3*i+2]);
			// Read extra data
			fread( temp, 1, this->extraVertex, fin ) ;
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
		// Read vertices in Ascii mode
		for ( i = 0 ; i < this->numVerts ; i ++ )
		{
			fgets ( temp, 1024, fin ) ;

			char seps[] = " ,\t\n";
			char* token = strtok( temp, seps ) ;

			sscanf( token, "%f", &(vertarray[ 3 * i ] ) ) ;
			token = strtok( NULL, seps ) ;
			sscanf( token, "%f", &(vertarray[ 3 * i + 1 ] ) ) ;
			token = strtok( NULL, seps ) ;
			sscanf( token, "%f", &(vertarray[ 3 * i + 2 ] ) ) ;
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
}


void PLYReader::writePLY(char *fname)
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
	// Write Vertices
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

void PLYReader::writeCrestLineInputTxt(char *fname)
{
	FILE *fout = fopen( fname, "w" );
	if ( fout == NULL )
	{
		printf("Unable to create file %s\n", fname);
		return;
	}

	// vertex properties
	fprintf( fout, "%d\n", numVerts);
	// face properties
	fprintf( fout, "%d\n", numTrians);
	// neighbor size
	int k = 1;
	fprintf( fout, "%d\n", k ) ;
	// flag for crest-line tracing (0 for no crest line tracing, 1 for yes)
	int crestlineflag = 1;
	fprintf( fout, "%d\n", crestlineflag ) ;

	////////////////////////////////////////////////////////////
	// Write Vertices
	float *v;
	for (int i = 0; i != numVerts; ++i)
	{
		v = getVertex(i);
		fprintf( fout, "%f %f %f\n",v[0], v[1], v[2]) ;
	}
	////////////////////////////////////////////////////////////
	// Write Faces
	int f[3];
	int *t;
	for ( int i = 0 ; i != numTrians; i++ )
	{
		t = getTriangle(i);
		for (int j = 0; j < 3; ++j)
		{
			f[j] = t[j];
		}
		fprintf( fout, "%d %d %d\n",f[0], f[1], f[2]);
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

void PLYReader::writePLY2(char *fname)
{
	FILE *fout = fopen( fname, "w" );
	if ( fout == NULL )
	{
		printf("Unable to create file %s\n", fname);
		return;
	}

	// vertex properties
	fprintf( fout, "%d\n", numVerts);
	// face properties
	fprintf( fout, "%d\n", numTrians);
	// Write Vertices
	float *v;
	for (int i = 0; i != numVerts; ++i)
	{
		v = getVertex(i);
		fprintf( fout, "%f %f %f\n",v[0], v[1], v[2]) ;
	}
	// Write Faces
	int f[3];
	int *t;
	for ( int i = 0 ; i != numTrians; i++ )
	{
		t = getTriangle(i);
		for (int j = 0; j < 3; ++j)
		{
			f[j] = t[j];
		}
		fprintf( fout, "3 ");
		fprintf( fout, "%d %d %d\n",f[0], f[1], f[2]);
	}
	fclose(fout);
}

void PLYReader::writeGTS(char *fname)
{
	/************************************************************************/
	/* 
	* Firstly get all edge pairs
	*/
	typedef map<pair<unsigned,unsigned>, unsigned> EdgePairType;
	EdgePairType edgePairs;
	EdgePairType::iterator iterSearch;
	vector<pair<unsigned,unsigned>> edgePairs_order;
	int *t;
	unsigned count = 1;
	for ( unsigned i = 0 ; i != numTrians; i++ )
	{
		t = getTriangle(i);
		for (int j = 0; j < 3; ++j)
		{
			unsigned a = t[j], b = t[(j+1)%3];
			if (a > b)
			{
				unsigned temp = b;
				b = a;
				a = temp;
			}
			iterSearch = edgePairs.find(pair<unsigned,unsigned>(a,b));
			if (iterSearch == edgePairs.end())
			{
				edgePairs[pair<unsigned,unsigned>(a,b)] = count;
				edgePairs_order.push_back(pair<unsigned,unsigned>(a,b));
				count++;
			}
		}
	}

	/************************************************************************/
	/* 
	* Next, write gts file
	*/

	FILE *fout = fopen( fname, "w" );
	if ( fout == NULL )
	{
		printf("Unable to create file %s\n", fname);
		return;
	}
	// Vertices, Edges and Faces number
	fprintf( fout, "%d %d %d\n", numVerts, edgePairs.size(), numTrians);
	
	/* 
	* Write vertices
	*/
	float *v;
	for (int i = 0; i != numVerts; ++i)
	{
		v = getVertex(i);
		fprintf( fout, "%f %f %f\n",v[0], v[1], v[2]) ;
	}
	/* 
	* Write edges
	*/
	for (vector<pair<unsigned,unsigned>>::iterator iter = edgePairs_order.begin(); 
		iter != edgePairs_order.end(); ++iter)
	{
		fprintf( fout, "%d %d\n", iter->first+1, iter->second+1);
	}
	/* 
	* Write faces
	*/
	for ( int i = 0 ; i != numTrians; i++ )
	{
		t = getTriangle(i);
		for (int j = 0; j < 3; ++j)
		{
			unsigned a = t[j], b = t[(j+1)%3];
			if (a > b)
			{
				unsigned temp = b;
				b = a;
				a = temp;
			}
			iterSearch = edgePairs.find(pair<unsigned,unsigned>(a,b));
			if (iterSearch != edgePairs.end())
			{
				fprintf( fout, "%d ", iterSearch->second);
			}
		}
		fprintf( fout, "\n");
	}
	fclose(fout);
	edgePairs.clear();
	edgePairs_order.clear();
}
