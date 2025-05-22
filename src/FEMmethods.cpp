#include "FEMmethods.h"

FiniteElementeBase::FiniteElementeBase(bool saveSolution) : saveSolution(saveSolution) {}

void FiniteElementeBase::updateElectrodeVoltage(int const& electrodeIndex,
    double const& voltage)
{

    // std::cout<<"setting electrode "<<electrodeIndex<<" to voltage "<<voltage<<"
    // indices: ";

    for (auto index : electrodeVertexIndices[electrodeIndex]) {
        (*solutionVector)[index] = voltage;
        // std::cout<<index<<" ";
    }
    // std::cout<<std::endl;
}

void FiniteElementeBase::initRun(bool initDevice /*= false*/)
{
    // see mfem ex1.cpp in mfem lib

    // 2.--> see see mfem ex1.cpp in mfem lib
    static bool mfem_initialized = false;
    initDevice = false;
    if (initDevice && !mfem_initialized) {
        Device device("cpu");
        mfem_initialized = true;
    }
    //    device.Print();

    // 6.--> see see mfem ex1.cpp in mfem lib

    if (mesh->bdr_attributes.Size()) {
        Array<int> ess_bdr(2);
        ess_bdr[0] = 0; // <<<<<<<<<<<<--------------  change this to 1 to make all
            // of boundaries essential
        ess_bdr[1] = 1;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    if (saveSolution) {
        std::ofstream mesh_ofs("finEle.mesh");
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);
    }
}

void FiniteElementeBase::run()
{

    // 7. Set up the linear form b(.) which corresponds to the right-hand side of
    //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
    //    the basis functions in the finite element fespace.
    b = new LinearForm(fespace);
    ConstantCoefficient one(1.0);
    ConstantCoefficient zero(0.0);
    b->AddDomainIntegrator(new DomainLFIntegrator(zero));
    b->Assemble();

    // 9. Set up the bilinear form a(.,.) on the finite element space
    //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
    //    domain integrator.
    a = new BilinearForm(fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));

    // 10. Assemble the bilinear form and the corresponding linear system,
    //     applying any necessary transformations such as: eliminating boundary
    //     conditions, applying conforming constraints for non-conforming AMR,
    //     static condensation, etc.
    a->Assemble();
    a->FormLinearSystem(ess_tdof_list, *solutionVector, *b, A, X, B);
    // 11.
    M = GSSmoother((SparseMatrix&)(*A));
    PCG(*A, M, B, X, 0, 1000, 1e-12, 0.0);

    // 12. Recover the solution as a finite element grid function.
    a->RecoverFEMSolution(X, *b, *solutionVector);

    // 13. Save the refined mesh and the solution. This output can be viewed later
    //    using GLVis: "glvis -m finEle.mesh -g laplace_solution.gf".
    if (saveSolution) {
        std::ofstream sol_ofs("laplace_solution" + std::to_string(runNumber) + ".gf");
        sol_ofs.precision(8);
        solutionVector->Save(sol_ofs);
        std::cout << "mfem run " << runNumber << " finished" << std::endl;
    }
    runNumber++;
    // *solutionVector = 0;

    delete a;
    delete b;
}

FiniteElementeBase::~FiniteElementeBase()
{
    delete fespace;
    delete fec;
    delete mesh;
}

FiniteElementeCircle::FiniteElementeCircle(double const& radius,
    int const& maxNumberOfElments,
    bool saveSolution)
    : FiniteElementeBase(saveSolution)
    , radius(radius)
{
    initMesh(maxNumberOfElments);

    // 5. --> see see mfem ex1.cpp in mfem lib
    fec = new H1_FECollection(order, dim);
    fespace = new FiniteElementSpace(mesh, fec);

    // 8.  --> see see mfem ex1.cpp in mfem lib
    solutionVector = new GridFunction(fespace); // changed to pointer, named x in example
    *solutionVector = 0;
}

void FiniteElementeCircle::initMesh(int const& maxNumberOfElments)
{

    layers = std::sqrt(maxNumberOfElments / 6.0);

    deltaR = radius / layers;

    int nv = 3 * (layers + 1) * layers + 1;
    int ne = 6 * layers * layers;
    int nb = 6 * layers;

    // std::cout<<"layers: "<<layers<<" elements "<<ne<<std::endl;

    mesh = new Mesh(2, nv, ne, nb, 2);

    // add vertices
    double R, phi;
    for (int l = layers; l > 0; l--) {
        R = l * deltaR;
        for (int i = 0; i < 6 * l; i++) {
            phi = 2 * PI * (i) / (6 * l);
            mesh->AddVertex(Vertex(R * cos(phi), R * sin(phi))());
        }
    }

    mesh->AddVertex(Vertex(double(0), double(0))());

    // add elements
    int triIndx[3];
    int innerIndex, outerIndex;
    for (int l = layers; l > 1; l--) { // last layer extra
        outerIndex = 0;
        innerIndex = 0;

        while (true) {
            triIndx[0] = innerIndex % (6 * (l - 1)) + (nv - 1) - 3 * l * (l - 1);
            triIndx[1] = outerIndex % (6 * l) + (nv - 1) - 3 * l * (l + 1);
            if ((outerIndex) / double(l) > (innerIndex) / double(l - 1)) {
                innerIndex++;
                triIndx[2] = innerIndex % (6 * (l - 1)) + (nv - 1) - 3 * l * (l - 1);
            } else {
                outerIndex++;
                triIndx[2] = outerIndex % (6 * l) + (nv - 1) - 3 * l * (l + 1);
            }
            // std::cout<<"out: "<<outerIndex<<" in: "<<innerIndex<<" tri:
            // "<<triIndx[0]<<" "<<triIndx[1]<<" "<<triIndx[2]<<std::endl;
            mesh->AddTri(triIndx);
            if (innerIndex % (6 * (l - 1)) == 0 and outerIndex % (6 * l) == 0) {
                break;
            }
        }
    }

    // add elements
    int tirIndx[3];
    for (int i = nv - 7; i < nv - 2; i++) {
        tirIndx[0] = i;
        tirIndx[1] = i + 1;
        tirIndx[2] = nv - 1;
        mesh->AddTri(tirIndx);
    }
    tirIndx[0] = nv - 2;
    tirIndx[1] = nv - 7;
    tirIndx[2] = nv - 1;
    mesh->AddTri(tirIndx);

    // add boundary
    int bdrIndx[2];
    for (int i = 0; i < 6 * layers - 1; i++) {
        bdrIndx[0] = i;
        bdrIndx[1] = i + 1;
        mesh->AddBdrSegment(bdrIndx);
    }
    bdrIndx[0] = 6 * layers - 1;
    bdrIndx[1] = 0;
    mesh->AddBdrSegment(bdrIndx);

    mesh->FinalizeTriMesh(0, 0, true);
}

void FiniteElementeCircle::setElectrode(double const& voltage, double begin,
    double end)
{

    // modulo  input values and round to nearest index
    int beginIdx = std::floor((begin / (2 * PI) - floor(begin / (2 * PI))) * 6 * layers);
    int endIdx = std::ceil((end / (2 * PI) - floor(end / (2 * PI))) * 6 * layers);

    // std::cout<<"beginIdx "<<beginIdx<< " endIdx "<<endIdx<<std::endl;

    // get vertex indices inside range
    if (beginIdx < endIdx) {
        electrodeVertexIndices.push_back(std::vector<int>(endIdx - beginIdx));
        std::iota(electrodeVertexIndices.back().begin(),
            electrodeVertexIndices.back().end(), beginIdx);
    } else {
        electrodeVertexIndices.push_back(
            std::vector<int>(6 * layers + endIdx - beginIdx));
        std::iota(electrodeVertexIndices.back().begin(),
            electrodeVertexIndices.back().begin() + 6 * layers - beginIdx,
            beginIdx);
        std::iota(electrodeVertexIndices.back().begin() + 6 * layers - beginIdx,
            electrodeVertexIndices.back().end(), 0);
    }

    // for(auto index:electrodeVertexIndices.back()){
    //     std::cout<<index<<" ";
    // }
    // std::cout<<std::endl;

    // set initial condition in solution vector
    for (auto index : electrodeVertexIndices.back()) {
        (*solutionVector)[index] = voltage;
    }

    // set boundary attribute to make (only) it mandatory
    for (auto i = electrodeVertexIndices.back().begin();
         i != electrodeVertexIndices.back().end() - 1; i++) {
        mesh->GetBdrElement(*i)->SetAttribute(2);
    }

    // move first and last vertex to exactly fit electrode boundaries
    if (electrodeVertexIndices.back().size() < 2) {
        throw std::length_error("could not find 2 vertices inside electrode range, "
                                "refine grid or enlarge electrode");
    } else {
        mesh->GetVertex(electrodeVertexIndices.back().front())[0] = radius * std::cos(begin);
        mesh->GetVertex(electrodeVertexIndices.back().front())[1] = radius * std::sin(begin);
        mesh->GetVertex(electrodeVertexIndices.back().back())[0] = radius * std::cos(end);
        mesh->GetVertex(electrodeVertexIndices.back().back())[1] = radius * std::sin(end);
    }

    numberOfElectrodes++;
}

double FiniteElementeCircle::getPotential(double const& x, double const& y)
{
    int layer = std::sqrt(x * x + y * y) / deltaR + 0.5;
    double phi = std::atan2(y, x) / (2 * PI);
    // std::cout<<"return index: "<<int(3*(layers+1)*layers-3*layer*(layer+1) +
    // (phi < 0 ? phi + 1 : phi) * 6 * layer + 0.5)<<std::endl;
    return (
        *solutionVector)[int(3 * (layers + 1) * layers - 3 * layer * (layer + 1) + (phi < 0 ? phi + 1 : phi) * 6 * layer + 0.5)];
}