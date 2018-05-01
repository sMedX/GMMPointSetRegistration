#include <itkEuler3DTransform.h>
#include <itkMesh.h>
#include <itkTransformMeshFilter.h>
#include <itkLBFGSOptimizer.h>

#include "itkInitializeTransform.h"
#include "itkInitializeRandomTransform.h"
#include "itkInitializeMetric.h"
#include "itkAddNoiseAndTransformMeshFilter.h"
#include "itkGMMPointSetToPointSetRegistrationMethod.h"
#include "itkPointSetPropertiesCalculator.h"
#include "itkPointSetToPointSetMetrics.h"

#include "itkIOutils.h"
#include "argsCustomParsers.h"

const unsigned int Dimension = 3;
typedef itk::Transform <double, Dimension, Dimension> TransformType;

typedef itk::Mesh<float, Dimension> MeshType;
typedef itk::PointSet<MeshType::PixelType, Dimension> PointSetType;

int main(int argc, char** argv) {

  // parse input arguments
  args::ArgumentParser parser("GMM-based point set to point set registration.", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  args::Group required(parser, "Required arguments:", args::Group::Validators::All);

  args::ValueFlag<std::string> argFixedFileName(required, "fixed", "The fixed mesh (point-set) file name", {'f', "fixed"});
  args::ValueFlag<std::string> argMovingFileName(required, "moving", "The moving mesh (point-set) file name", {'m', "moving"});

  args::ValueFlag<size_t> argNumberOfLevels(parser, "levels", "The number of levels", {"levels"});
  args::ValueFlag<size_t> argNumberOfIterations(parser, "iterations", "The number of iterations", {"iterations"}, 1000);
  args::ValueFlag<size_t> argNumberOfEvaluations(parser, "evaluations", "The number of evaluations", {"evaluations"}, 1000);

  args::Flag trace(parser, "trace", "Optimizer iterations tracing", {"trace"});

  const std::string transformDescription =
    "The type of transform (That is number):\n"
    "  0 : Translation\n"
    "  1 : Versor3D\n"
    "  2 : Similarity\n"
    "  3 : ScaleSkewVersor3D\n";

  args::ValueFlag<size_t> argTypeOfTransform(parser, "transform", transformDescription, {"transform"}, 0);

  const std::string metricDescription =
    "The type of metric (That is number):\n"
    "  0 : L2Rigid\n"
    "  1 : L2\n"
    "  2 : KC\n";

  args::ValueFlag<size_t> argTypeOfMetric(parser, "metric", metricDescription, {"metric"}, 0);

  try {
    parser.ParseCLI(argc, argv);
  }
  catch (args::Help) {
    std::cout << parser;
    return EXIT_SUCCESS;
  }
  catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return EXIT_FAILURE;
  }
  catch (args::ValidationError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return EXIT_FAILURE;
  }

  std::string fixedFileName = args::get(argFixedFileName);
  std::string movingFileName = args::get(argMovingFileName);
  size_t numberOfEvaluations = args::get(argNumberOfEvaluations);
  size_t numberOfLevels = args::get(argNumberOfLevels);
  size_t numberOfIterations = args::get(argNumberOfIterations);
  size_t typeOfTransform = args::get(argTypeOfTransform);
  size_t typeOfMetric = args::get(argTypeOfMetric);

  std::cout << "options" << std::endl;
  std::cout << "number of evaluations " << argNumberOfEvaluations << std::endl;
  std::cout << "number of levels      " << numberOfLevels << std::endl;
  std::cout << "number of iterations  " << numberOfIterations << std::endl;
  std::cout << "transform " << typeOfTransform << std::endl;
  std::cout << "metric    " << typeOfMetric << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  // read meshes
  MeshType::Pointer fixedMesh = MeshType::New();
  if (!readMesh<MeshType>(fixedMesh, fixedFileName)) {
    return EXIT_FAILURE;
  }

  std::cout << "fixed mesh " << fixedFileName << std::endl;
  std::cout << "number of points " << fixedMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  MeshType::Pointer initialMovingMesh = MeshType::New();
  if (!readMesh<MeshType>(initialMovingMesh, movingFileName)) {
    return EXIT_FAILURE;
  }

  std::cout << "moving mesh " << movingFileName << std::endl;
  std::cout << "number of points " << initialMovingMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  // random transform initializer
  typedef itk::PointSetPropertiesCalculator<PointSetType> MovingPointSetPropertiesCalculatorType;
  MovingPointSetPropertiesCalculatorType::Pointer movingPointSetCalculator = MovingPointSetPropertiesCalculatorType::New();
  movingPointSetCalculator->SetPointSet(initialMovingMesh);
  try 
  {
    movingPointSetCalculator->Compute();
  }
  catch (itk::ExceptionObject& excep) 
  {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  movingPointSetCalculator->Print(std::cout);

  typedef itk::InitializeRandomTransform<double> RandomTransformInitializerType;
  RandomTransformInitializerType::Pointer initializerRandomTransform = RandomTransformInitializerType::New();
  initializerRandomTransform->SetCenter(movingPointSetCalculator->GetCenter());
  initializerRandomTransform->SetTypeOfTransform(typeOfTransform);

  // initialize metric
  typedef itk::InitializeMetric<PointSetType, PointSetType> InitializeMetricType;
  InitializeMetricType::Pointer initializerMetric = InitializeMetricType::New();
  initializerMetric->SetTypeOfMetric(typeOfMetric);
  try {
    initializerMetric->Initialize();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  initializerMetric->Print(std::cout);
  InitializeMetricType::MetricType::Pointer metric = initializerMetric->GetMetric();

  // initialize optimizer
  typedef itk::LBFGSOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetMaximumNumberOfFunctionEvaluations(numberOfIterations);
  optimizer->SetTrace(trace);

  //--------------------------------------------------------------------
  // compute statistics
  for (size_t count = 0; count < numberOfEvaluations; ++count) 
  {
    std::cout << "---------------------------------------------------------------------" << std::endl;
    std::cout << "Perform registration: " << count+1 << " (" << numberOfEvaluations << ")" << std::endl;
    std::cout << std::endl;

    // generate random transform 
    try {
      initializerRandomTransform->Initialize();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    initializerRandomTransform->Print(std::cout);

    // transform initial moving mesh
    typedef itk::AddNoiseAndTransformMeshFilter<MeshType, MeshType, TransformType> RandomTransformFilterType;
    RandomTransformFilterType::Pointer transformMesh1 = RandomTransformFilterType::New();
    transformMesh1->SetInput(initialMovingMesh);
    transformMesh1->SetTransform(initializerRandomTransform->GetTransform());
    transformMesh1->SetStandardDeviation(1);
    try {
      transformMesh1->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    MeshType::Pointer movingMesh = transformMesh1->GetOutput();

    // initialize initial registration
    typedef itk::PointSetPropertiesCalculator<PointSetType> MovingPointSetPropertiesCalculatorType;
    MovingPointSetPropertiesCalculatorType::Pointer movingPointSetCalculator = MovingPointSetPropertiesCalculatorType::New();
    movingPointSetCalculator->SetPointSet(movingMesh);
    try {
      movingPointSetCalculator->Compute();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    movingPointSetCalculator->Print(std::cout);

    typedef itk::InitializeTransform<double> TransformInitializerType;
    TransformInitializerType::Pointer initializerTransform = TransformInitializerType::New();
    initializerTransform->SetFixedLandmark(movingPointSetCalculator->GetCenter());
    initializerTransform->SetMovingLandmark(movingPointSetCalculator->GetCenter());
    initializerTransform->SetTypeOfTransform(typeOfTransform);
    try {
      initializerTransform->Initialize();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    initializerTransform->Print(std::cout);
    TransformType::Pointer transform = initializerTransform->GetTransform();

    optimizer->SetScales(initializerTransform->GetScales());

    // perform registration
    typedef itk::GMMPointSetToPointSetRegistrationMethod<PointSetType> GMMPointSetToPointSetRegistrationMethodType;
    GMMPointSetToPointSetRegistrationMethodType::Pointer registration = GMMPointSetToPointSetRegistrationMethodType::New();
    registration->SetFixedPointSet(fixedMesh->GetPoints());
    registration->SetMovingPointSet(transformMesh1->GetOutput());
    registration->SetMetric(initializerMetric->GetMetric());
    registration->SetTransform(transform);
    registration->SetNumberOfLevels(numberOfLevels);
    registration->SetOptimizer(optimizer);
    try {
      registration->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }

    // transform moving mesh
    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformMeshFilterType;
    TransformMeshFilterType::Pointer transformMesh2 = TransformMeshFilterType::New();
    transformMesh2->SetInput(movingMesh);
    transformMesh2->SetTransform(transform);
    try {
      transformMesh2->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }

    // compute metrics
    typedef itk::PointSetToPointSetMetrics<PointSetType> PointSetToPointSetMetricsType;
    PointSetToPointSetMetricsType::Pointer metrics = PointSetToPointSetMetricsType::New();
    metrics->SetFixedPointSet(fixedMesh->GetPoints());
    metrics->SetTargetPointSet(initialMovingMesh->GetPoints());
    metrics->SetMovingPointSet(transformMesh1->GetOutput());
    metrics->Compute();
    metrics->PrintReport(std::cout);

    metrics->SetFixedPointSet(fixedMesh->GetPoints());
    metrics->SetTargetPointSet(initialMovingMesh->GetPoints());
    metrics->SetMovingPointSet(transformMesh2->GetOutput());
    metrics->Compute();
    metrics->PrintReport(std::cout);
  }

  return EXIT_SUCCESS;
}
