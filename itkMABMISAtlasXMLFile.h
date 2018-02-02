#ifndef __itkMABMISAtlasXMLFile_h
#define __itkMABMISAtlasXMLFile_h

#include "itkXMLFile.h"
#include <string>
#include <vector>

#ifdef _WIN32
static const char FILESEP = '\\';
#else
static const char FILESEP = '/';
#endif

namespace itk
{
class MABMISImageData
{
public:
  std::string              m_DataDirectory;
  std::string              m_OutputDirectory;
  int                      m_NumberImageData;
  std::vector<std::string> m_ImageFileNames;
  std::vector<std::string> m_SegmentationFileNames;

  MABMISImageData():
    m_DataDirectory(""),
    m_OutputDirectory(""),
    m_NumberImageData(0),
    m_ImageFileNames(0),
    m_SegmentationFileNames(0)
  {
  }
  ~MABMISImageData()
  {
  }
};

class MABMISAtlas
{
public:
  std::string m_AtlasDirectory;

  int                      m_NumberAllAtlases;
  int                      m_NumberSimulatedAtlases;
  std::vector<std::string> m_AtlasFilenames;
  std::vector<std::string> m_AtlasSegmentationFilenames;

  std::vector<bool> m_IsSimulatedImage;

  std::vector<int> m_RealImageIDs;
  std::vector<int> m_SimulatedImageIDs;

  std::vector<int> m_Tree;
  int              m_TreeSize;
  int              m_TreeRoot;
  int              m_TreeHeight;

  MABMISAtlas()
  {
    m_NumberAllAtlases = 0;
    m_NumberSimulatedAtlases = 0;
    m_AtlasFilenames.resize(0);

    m_Tree.resize(0);
    m_TreeSize = 0;
    m_TreeRoot = -1;
    m_TreeHeight = 0;
  }
  ~MABMISAtlas()
  {
  }
};

class MABMISImageDataXMLFileReader :
  public XMLReader<MABMISImageData>
{
public:
  /** Standard typedefs */
  typedef MABMISImageDataXMLFileReader Self;
  typedef XMLReader<MABMISImageData>   Superclass;
  typedef SmartPointer<Self>           Pointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISImageDataXMLFileReader, XMLReader);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
public:
  /** Determine if a file can be read */
  virtual int CanReadFile(const char *name);

protected:
  MABMISImageDataXMLFileReader()
  {
    m_OutputObject = new itk::MABMISImageData;
    m_ImageCount = 0;
    m_ImageData = 0;
  }
  virtual ~MABMISImageDataXMLFileReader()
  {
  }

  virtual void StartElement(const char *name, const char * *atts);

  virtual void EndElement(const char *name);

  virtual void CharacterDataHandler(const char *inData, int inLength);

private:
  MABMISImageDataXMLFileReader(const Self &) ITK_DELETED_FUNCTION;
  void operator=(const Self &) ITK_DELETED_FUNCTION;
  MABMISImageData* m_ImageData;
  std::string      m_CurCharacterData;

  int m_ImageCount;
};

class MABMISImageDataXMLFileWriter :
  public XMLWriterBase<MABMISImageData>
{
public:
  /** standard typedefs */
  typedef XMLWriterBase<MABMISImageData> Superclass;
  typedef MABMISImageDataXMLFileWriter   Self;
  typedef SmartPointer<Self>             Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISImageDataXMLFileWriter,  XMLWriterBase<MABMISImageData> );

  /** Test whether a file is writable. */
  virtual int CanWriteFile(const char *name);

  /** Actually write out the file in question */
  int WriteFile() ITK_OVERRIDE;

protected:
  MABMISImageDataXMLFileWriter()
  {
  }

  virtual ~MABMISImageDataXMLFileWriter()
  {
  }

private:
  MABMISImageDataXMLFileWriter(const Self &) ITK_DELETED_FUNCTION;
  void operator=(const Self &) ITK_DELETED_FUNCTION;
};

class MABMISAtlasXMLFileReader :
  public XMLReader<MABMISAtlas>
{
public:
  /** Standard typedefs */
  typedef MABMISAtlasXMLFileReader      Self;
  typedef XMLReader<MABMISImageData>    Superclass;
  typedef SmartPointer<Self>            Pointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISAtlasXMLFileReader, XMLReader);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  virtual void GenerateOutputInformation();

public:
  /** Determine if a file can be read */
  virtual int CanReadFile(const char *name);

protected:
  MABMISAtlasXMLFileReader()
  {
    m_OutputObject = new itk::MABMISAtlas;
    m_Atlas = 0;
  }
  virtual ~MABMISAtlasXMLFileReader()
  {
  }

  virtual void StartElement(const char *name, const char * *atts);

  virtual void EndElement(const char *name);

  virtual void CharacterDataHandler(const char *inData, int inLength);

private:
  MABMISAtlasXMLFileReader(const Self &) ITK_DELETED_FUNCTION;
  void operator=(const Self &) ITK_DELETED_FUNCTION;
  MABMISAtlas* m_Atlas;
  std::string  m_CurCharacterData;
};

class MABMISAtlasXMLFileWriter :
  public XMLWriterBase<MABMISAtlas>
{
public:
  /** standard typedefs */
  typedef XMLWriterBase<MABMISAtlas> Superclass;
  typedef MABMISAtlasXMLFileWriter   Self;
  typedef SmartPointer<Self>         Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISAtlasXMLFileWriter,  XMLWriterBase<MABMISAtlas> );

  /** Test whether a file is writable. */
  virtual int CanWriteFile(const char *name);

  /** Actually write out the file in question */
  int WriteFile() ITK_OVERRIDE;

protected:
  MABMISAtlasXMLFileWriter()
  {
  }

  virtual ~MABMISAtlasXMLFileWriter()
  {
  }

private:
  MABMISAtlasXMLFileWriter(const Self &) ITK_DELETED_FUNCTION;
  void operator=(const Self &) ITK_DELETED_FUNCTION;
};
}

#endif
