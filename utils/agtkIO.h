#ifndef __agtkIO_h
#define __agtkIO_h

#include <iostream>

#include <itkMeshFileReader.h>
#include <itkMeshFileWriter.h>


namespace agtk
{

//! Reads a mesh from a file
template <typename TMesh>
bool readMesh(typename TMesh::Pointer mesh, const std::string& fileName)
{
  typedef itk::MeshFileReader<TMesh> MeshFileReader;
  typename MeshFileReader::Pointer reader = MeshFileReader::New();
  reader->SetFileName(fileName);

  try {
    reader->Update();
  }
  catch (itk::ExceptionObject& err) {
    std::cerr << "Unable to read mesh from file '" << fileName << "'" << std::endl;
    std::cerr << "Error: " << err << std::endl;
    return false;
  }

  mesh->Graft(reader->GetOutput());
  return true;
}

//! Writes a mesh to a file
template <typename TMesh>
bool writeMesh(const TMesh* mesh, const std::string& fileName)
{
  typedef itk::MeshFileWriter<TMesh> MeshFileWriter;
  typename MeshFileWriter::Pointer writer = MeshFileWriter::New();

  writer->SetFileName(fileName);
  writer->SetInput(mesh);
  writer->SetFileTypeAsASCII();
  try {
    writer->Update();
  }
  catch (itk::ExceptionObject& err) {
    std::cerr << "Unable to write mesh to file '" << fileName << "'" << std::endl;
    std::cerr << "Error: " << err << std::endl;
    return false;
  }

  return true;
}
}

#endif // __agtkIO_h
