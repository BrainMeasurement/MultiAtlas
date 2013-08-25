#ifndef __itkMABMISAtlasXMLFile_cxx
#define __itkMABMISAtlasXMLFile_cxx

#include "itksys/SystemTools.hxx"
#include "itkMABMISAtlasXMLFile.h"

#define RAISE_EXCEPTION(s)                                     \
    { \
    ExceptionObject \
    exception( \
      __FILE__, \
      __LINE__); \
    exception. \
    SetDescription( \
      s);                     \
    throw exception; \
    }

namespace itk
{
///------------------------------------------------------------
// ----  class MABMISImageDataXMLFileReader  ----------------------
///------------------------------------------------------------

int
MABMISImageDataXMLFileReader::CanReadFile(const char *name)
{
  if( !itksys::SystemTools::FileExists(name)
      || itksys::SystemTools::FileIsDirectory(name)
      || itksys::SystemTools::FileLength(name) == 0 )
    {
    return 0;
    }
  return 1;
}

void
MABMISImageDataXMLFileReader::StartElement( const char *name, const char * *atts )
{
  if( itksys::SystemTools::Strucmp(name, "MABMISImageData") == 0 )
    {
    m_ImageData = new itk::MABMISImageData;
    }
  else if( itksys::SystemTools::Strucmp(name, "DataPath") == 0 )
    {
    }
  else if( itksys::SystemTools::Strucmp(name, "NumberOfImageData") == 0 )
    {
    }
  else if( itksys::SystemTools::Strucmp(name, "Dataset") == 0 )
    {
    int         i = 0;
    std::string imageName("");
    std::string segImageName("");
    while( atts[i] )
      {
      const std::string key(atts[i]);
      ++i;
      const std::string value(atts[i]);
      ++i;
      if( itksys::SystemTools::Strucmp(key.c_str(), "ImageFile") == 0 )
        {
        imageName = value;
        }
      else if( itksys::SystemTools::Strucmp(key.c_str(), "SegmentationFile") == 0 )
        {
        segImageName = value;
        }
      }

    if( this->m_ImageCount < m_ImageData->m_NumberImageData )
      {
      m_ImageData->m_ImageFileNames[this->m_ImageCount] = imageName;
      m_ImageData->m_SegmentationFileNames[this->m_ImageCount] = segImageName;

      this->m_ImageCount++;
      }
    }

  return;
}

void
MABMISImageDataXMLFileReader::EndElement(const char *name)
{
  if( itksys::SystemTools::Strucmp(name, "MABMISImageData") == 0 )
    {
    delete m_OutputObject;
    m_OutputObject = &(*(this->m_ImageData ) );
    }
  else
    {
    if( m_ImageData == 0 )
      {
      RAISE_EXCEPTION("The xml format is not right. Please check it!");
      }
    }
  if( itksys::SystemTools::Strucmp(name, "DataPath") == 0 )
    {
    this->m_ImageData->m_DataDirectory = this->m_CurCharacterData;
    }
  if( itksys::SystemTools::Strucmp(name, "NumberOfImageData") == 0 )
    {
    int numData = atoi(this->m_CurCharacterData.c_str() );
    this->m_ImageData->m_NumberImageData = numData;
    this->m_ImageData->m_ImageFileNames.resize(numData);
    this->m_ImageData->m_SegmentationFileNames.resize(numData);
    this->m_ImageCount = 0;
    }
}

static std::string StripLastNewline(const std::string & input)
{
  std::string output = input;
  size_t      last_element = input.size() - 1;

  if( input[last_element] == '\n' )
    {
    output = input.substr(0, last_element);
    }
  return output;
}

void
MABMISImageDataXMLFileReader::CharacterDataHandler(const char *inData, int inLength)
{
  this->m_CurCharacterData = "";
  for( int i = 0; i < inLength; ++i )
    {
    m_CurCharacterData = m_CurCharacterData + inData[i];
    }
  m_CurCharacterData = StripLastNewline(m_CurCharacterData);
}

///------------------------------------------------------------
// class MABMISAtlasXMLFileReader  ----------------------------
///------------------------------------------------------------

int
MABMISAtlasXMLFileReader::CanReadFile(const char *name)
{
  if( !itksys::SystemTools::FileExists(name)
      || itksys::SystemTools::FileIsDirectory(name)
      || itksys::SystemTools::FileLength(name) == 0 )
    {
    return 0;
    }
  return 1;
}

void
MABMISAtlasXMLFileReader::GenerateOutputInformation()
{
  this->parse();

  // validate the results are right.
  int numRealImages = this->m_OutputObject->m_NumberAllAtlases - this->m_OutputObject->m_NumberSimulatedAtlases;
  int numSimImages = this->m_OutputObject->m_NumberSimulatedAtlases;
  this->m_OutputObject->m_SimulatedImageIDs.resize(numSimImages);
  this->m_OutputObject->m_RealImageIDs.resize(numRealImages);

  int countRealImages = 0, countSimImages = 0;
  for( int n = 0; n < this->m_OutputObject->m_IsSimulatedImage.size(); n++ )
    {
    if( this->m_OutputObject->m_IsSimulatedImage[n] )
      {
      if( countSimImages >= numSimImages ) // means info in xml file is inconsistent
        {
        continue;           // need do something: later
        }
      this->m_OutputObject->m_SimulatedImageIDs[countSimImages] = n;
      countSimImages++;
      }
    else
      {
      if( countRealImages >= numRealImages ) // means info in xml file is inconsistent
        {
        continue;           // need do something
        }
      this->m_OutputObject->m_RealImageIDs[countRealImages] = n;
      countRealImages++;
      }
    }
}

void
MABMISAtlasXMLFileReader::StartElement( const char *name, const char * *atts )
{
  if( itksys::SystemTools::Strucmp(name, "MABMISATLAS") == 0 )
    {
    m_Atlas = new itk::MABMISAtlas;
    }
  else if( itksys::SystemTools::Strucmp(name, "ATLASFILELIST") == 0 )
    {
    m_Atlas->m_AtlasFilenames.resize(m_Atlas->m_NumberAllAtlases);
    m_Atlas->m_AtlasSegmentationFilenames.resize(m_Atlas->m_NumberAllAtlases);
    m_Atlas->m_IsSimulatedImage.resize(m_Atlas->m_NumberAllAtlases);
    }
  else if( itksys::SystemTools::Strucmp(name, "ATLAS") == 0 )
    {
    int         i = 0;
    int         id;
    std::string fname = "";
    std::string segFName = "";
    bool        isSimulated = false;

    while( atts[i] )
      {
      std::string key(atts[i]);
      ++i;
      std::string value(atts[i]);
      ++i;

      if( itksys::SystemTools::Strucmp(key.c_str(), "ID") == 0 )
        {
        id = atoi(value.c_str() );
        }
      else if( itksys::SystemTools::Strucmp(key.c_str(), "FILENAME") == 0 )
        {
        fname = value;
        }
      else if( itksys::SystemTools::Strucmp(key.c_str(), "SEGMENTATIONFILENAME") == 0 )
        {
        segFName = value;
        }
      else if( itksys::SystemTools::Strucmp(key.c_str(), "SIMULATED") == 0 )
        {
        if( itksys::SystemTools::Strucmp(value.c_str(), "TRUE") ||
            itksys::SystemTools::Strucmp(value.c_str(), "1") ||
            itksys::SystemTools::Strucmp(value.c_str(), "T") )
          {
          isSimulated = true;
          }
        else
          {
          isSimulated = false;
          }
        }
      }

    m_Atlas->m_AtlasFilenames[id] = fname;
    m_Atlas->m_AtlasSegmentationFilenames[id] = segFName;
    m_Atlas->m_IsSimulatedImage[id] = isSimulated;
    }
  else if( itksys::SystemTools::Strucmp(name, "ATLASFILELIST") == 0 )
    {
    m_Atlas->m_AtlasFilenames.resize(m_Atlas->m_NumberAllAtlases);
    }
  else if( itksys::SystemTools::Strucmp(name, "Tree") == 0 )
    {
    m_Atlas->m_Tree.resize(m_Atlas->m_NumberAllAtlases);
    for( int n = 0; n < m_Atlas->m_Tree.size(); n++ )
      {
      m_Atlas->m_Tree[n] = -1;
      }
    }
  else if( itksys::SystemTools::Strucmp(name, "Node") == 0 )
    {
    int i = 0;
    int id, parentID;
    ;
    while( atts[i] )
      {
      const std::string key(atts[i]);
      ++i;
      const std::string value(atts[i]);
      ++i;
      if( itksys::SystemTools::Strucmp(key.c_str(), "ID") == 0 )
        {
        id = atoi(value.c_str() );
        }
      else if( itksys::SystemTools::Strucmp(key.c_str(), "ParentID") == 0 )
        {
        parentID = atoi(value.c_str() );
        }
      }

    m_Atlas->m_Tree[id] = parentID;
    }
}

void
MABMISAtlasXMLFileReader::EndElement(const char *name)
{
  if( itksys::SystemTools::Strucmp(name, "MABMISAtlas") == 0 )
    {
    delete m_OutputObject;
    m_OutputObject = &(*(this->m_Atlas ) );
    }
  else
    {
    if( m_Atlas == 0 )
      {
      RAISE_EXCEPTION("The xml format is not right");
      }
    }
  if( itksys::SystemTools::Strucmp(name, "AtlasDirectory") == 0 )
    {
    this->m_Atlas->m_AtlasDirectory = this->m_CurCharacterData;
    }
  if( itksys::SystemTools::Strucmp(name, "NumberOfAtlases") == 0 )
    {
    int numAtlas = atoi(this->m_CurCharacterData.c_str() );
    this->m_Atlas->m_NumberAllAtlases = numAtlas;
    }
  if( itksys::SystemTools::Strucmp(name, "NumberOfSimulatedAtlases") == 0 )
    {
    int numSimAtlas = atoi(this->m_CurCharacterData.c_str() );
    this->m_Atlas->m_NumberSimulatedAtlases = numSimAtlas;
    }
  if( itksys::SystemTools::Strucmp(name, "TreeSize") == 0 )
    {
    int treeSize = atoi(this->m_CurCharacterData.c_str() );
    this->m_Atlas->m_TreeSize = treeSize;
    }
  // if ( itksys::SystemTools::Strucmp(name, "TreeHeight") == 0 )
  // {
  //	int treeHeight = atoi(this->m_CurCharacterData.c_str());
  //	this->m_Atlas->m_TreeHeight = treeHeight;
  // }
  // if ( itksys::SystemTools::Strucmp(name, "TreeRootID") == 0 )
  // {
  //	int root = atoi(this->m_CurCharacterData.c_str());
  //	this->m_Atlas->m_TreeRoot = root;
  // }
}

void
MABMISAtlasXMLFileReader::CharacterDataHandler(const char *inData, int inLength)
{
  this->m_CurCharacterData = "";
  for( int i = 0; i < inLength; ++i )
    {
    m_CurCharacterData = m_CurCharacterData + inData[i];
    }
  m_CurCharacterData = StripLastNewline(m_CurCharacterData);
}

///------------------------------------------------------------
// ----  class MABMISAtlasXMLFileWriter  ----------------------
///------------------------------------------------------------

int
MABMISAtlasXMLFileWriter::CanWriteFile( const char *itkNotUsed(name) )
{
  return true;                  // not sure what else to say
}

int
MABMISAtlasXMLFileWriter::WriteFile()
{
  // sanity checks
  if( m_InputObject == 0 )
    {
    std::string errmsg("No MABMISAtlas to Write");
    RAISE_EXCEPTION(errmsg);
    }
  if( m_Filename.length() == 0 )
    {
    std::string errmsg("No filename given");
    RAISE_EXCEPTION(errmsg);
    }
  std::ofstream output( m_Filename.c_str() );
  if( output.fail() )
    {
    std::string errmsg("Can't Open ");
    errmsg += m_Filename;
    RAISE_EXCEPTION(errmsg);
    }

  WriteStartElement("?xml version=\"1.0\"?", output);
  output << std::endl;
  WriteStartElement("!DOCTYPE MABMISAtlas", output);
  output << std::endl;

  // Write out metadata

  WriteStartElement("MABMISAtlas", output);
  output << std::endl;

  // Path to atlas files
  // leave as blank for now
  WriteStartElement("AtlasDirectory", output);
  output << this->m_InputObject->m_AtlasDirectory;
  WriteEndElement("AtlasDirectory", output);
  output << std::endl;

  // number of atlas
  WriteStartElement("NumberOfAtlases", output);
  output << this->m_InputObject->m_NumberAllAtlases;
  WriteEndElement("NumberOfAtlases", output);
  output << std::endl;
  // number of atlas
  WriteStartElement("NumberOfSimulatedAtlases", output);
  output << this->m_InputObject->m_NumberSimulatedAtlases;
  WriteEndElement("NumberOfSimulatedAtlases", output);
  output << std::endl;

  // atlas file list
  WriteStartElement("AtlasFileList", output);
  output << std::endl;
  for( int n = 0; n < this->m_InputObject->m_AtlasFilenames.size(); n++ )
    {
    output << "<Atlas ";
    output << "ID=\"" << n << "\"";
    output << " FileName=\"" << this->m_InputObject->m_AtlasFilenames[n] << "\"";
    if( !this->m_InputObject->m_IsSimulatedImage[n] )
      {
      output << " SegmentationFileName=\"";
      output << this->m_InputObject->m_AtlasSegmentationFilenames[n] << "\"";
      output << " Simulated=\"false\"";
      }
    else
      {
      output << " Simulated=\"true\"";
      }

    output << "/>" << std::endl;
    }

  WriteEndElement("AtlasFileList", output);
  output << std::endl;

  // tree
  WriteStartElement("AtlasTree", output);
  output << std::endl;
  // tree size
  WriteStartElement("TreeSize", output);
  output << this->m_InputObject->m_Tree.size();
  WriteEndElement("TreeSize", output);
  output << std::endl;
  // tree height
  // WriteStartElement("TreeHeight", output);
  // output << this->m_InputObject->m_TreeHeight;
  // WriteEndElement("TreeHeight", output);
  // output<< std::endl;

  // tree root ID
  // WriteStartElement("TreeRootID", output);
  // output << this->m_InputObject->m_TreeRoot;
  // WriteEndElement("TreeRootID", output);
  // output<< std::endl;
  // write the tree
  WriteStartElement("Tree", output);
  output << std::endl;
  for( int i = 0; i < this->m_InputObject->m_Tree.size(); ++i )
    {
    output << "<Node ";
    output << "ID=\"" << i << "\"";
    output << " ParentID=\"" << this->m_InputObject->m_Tree[i] << "\"";
    output << "/>" << std::endl;
    }
  WriteEndElement("Tree", output);
  output << std::endl;

  WriteEndElement("AtlasTree", output);
  output << std::endl;

  WriteEndElement("MABMISAtlas", output);
  output << std::endl;
  output.close();

  return 0;
}
}

#endif
