# Geant4: From Particle Collisions to Molecular Interactions - API Implementation Guide

## Table of Contents
1. [Introduction and Core APIs](#introduction-and-core-apis)
2. [From Accelerator Physics to Atomic Scale - Implementation](#from-accelerator-physics-to-atomic-scale---implementation)
3. [Atomic Particles to Molecular Systems - APIs](#atomic-particles-to-molecular-systems---apis)
4. [Interaction Physics Implementation](#interaction-physics-implementation)
5. [Physics Lists and Model APIs](#physics-lists-and-model-apis)
6. [Bindings and Extensions](#bindings-and-extensions)
7. [Advanced API Usage and Custom Physics](#advanced-api-usage-and-custom-physics)

## Introduction and Core APIs

**Geant4** (GEometry ANd Tracking) is a comprehensive Monte Carlo simulation toolkit with extensive C++ APIs for modeling particle-matter interactions across scales from TeV to eV.

### Core Architecture Classes

```cpp
// Essential Geant4 manager classes
#include "G4RunManager.hh"           // Main simulation manager
#include "G4MTRunManager.hh"          // Multi-threaded version
#include "G4VisExecutive.hh"          // Visualization manager
#include "G4UImanager.hh"             // User interface manager
#include "G4SteppingManager.hh"       // Controls particle stepping
#include "G4EventManager.hh"          // Event processing
#include "G4TrackingManager.hh"       // Track management
```

### Mandatory User Classes

Every Geant4 application must implement three abstract classes:

```cpp
// 1. Detector Construction
class MyDetector : public G4VUserDetectorConstruction {
public:
    virtual G4VPhysicalVolume* Construct() override;
    virtual void ConstructSDandField() override;  // Sensitive detectors
};

// 2. Physics List
class MyPhysics : public G4VModularPhysicsList {
public:
    MyPhysics();
    virtual void SetCuts() override;
};

// 3. Primary Generator
class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
    virtual void GeneratePrimaries(G4Event* event) override;
};
```

### Initialization Pattern

```cpp
int main() {
    // Choose appropriate RunManager
    #ifdef G4MULTITHREADED
        G4MTRunManager* runManager = new G4MTRunManager;
        runManager->SetNumberOfThreads(8);
    #else
        G4RunManager* runManager = new G4RunManager;
    #endif
    
    // Register mandatory classes
    runManager->SetUserInitialization(new MyDetector());
    runManager->SetUserInitialization(new MyPhysics());
    runManager->SetUserAction(new MyPrimaryGenerator());
    
    // Initialize kernel
    runManager->Initialize();
    
    // Start simulation
    runManager->BeamOn(numberOfEvents);
}
```

## From Accelerator Physics to Atomic Scale - Implementation

### High-Energy Collisions APIs (TeV → GeV)

#### Particle Definition APIs

```cpp
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"

// Accessing particle properties
G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
G4ParticleDefinition* proton = particleTable->FindParticle("proton");
G4ParticleDefinition* electron = particleTable->FindParticle("e-");

// Creating exotic particles
#include "G4ParticlePropertyTable.hh"
G4ParticleDefinition* higgs = particleTable->FindParticle(25); // PDG code

// Dynamic particles with momentum
G4DynamicParticle* dynamicProton = new G4DynamicParticle(
    proton,
    G4ThreeVector(0, 0, 1),  // Direction
    7*TeV                      // Kinetic energy
);
```

#### String and QCD Model Implementation

```cpp
#include "G4FTFModel.hh"              // Fritiof string model
#include "G4ExcitedStringDecay.hh"   // String fragmentation
#include "G4QGSMFragmentation.hh"     // Quark-gluon string fragmentation
#include "G4LundStringFragmentation.hh"

// Setting up high-energy hadronic models
class HighEnergyPhysics : public G4VPhysicsConstructor {
public:
    void ConstructProcess() override {
        // Get FTF (Fritiof) model
        G4FTFModel* ftf = new G4FTFModel();
        
        // Configure string fragmentation
        G4ExcitedStringDecay* stringDecay = 
            new G4ExcitedStringDecay(new G4LundStringFragmentation());
        ftf->SetFragmentationModel(stringDecay);
        
        // Set high-energy generator
        G4TheoFSGenerator* theo = new G4TheoFSGenerator("FTFP");
        theo->SetHighEnergyGenerator(ftf);
        
        // Apply to protons above 4 GeV
        G4ProcessManager* pManager = G4Proton::Proton()->GetProcessManager();
        G4HadronInelasticProcess* inelProcess = new G4HadronInelasticProcess();
        inelProcess->RegisterMe(theo);
        pManager->AddDiscreteProcess(inelProcess);
    }
};
```

### Intermediate Energy Cascade APIs (GeV → MeV)

#### Bertini Cascade Implementation

```cpp
#include "G4CascadeInterface.hh"      // Bertini cascade
#include "G4BinaryCascade.hh"          // Binary cascade
#include "G4PreCompoundModel.hh"       // Pre-equilibrium

// Configure cascade models
void SetupCascadePhysics() {
    // Bertini cascade (preferred for 0-10 GeV)
    G4CascadeInterface* bertini = new G4CascadeInterface();
    bertini->SetMinEnergy(0.0*GeV);
    bertini->SetMaxEnergy(10.0*GeV);
    
    // Binary cascade (alternative, better for ions)
    G4BinaryCascade* binary = new G4BinaryCascade();
    binary->SetMinEnergy(0.0*GeV);
    binary->SetMaxEnergy(1.5*GeV);
    
    // Pre-compound/de-excitation
    G4PreCompoundModel* preco = new G4PreCompoundModel();
    G4ExcitationHandler* handler = new G4ExcitationHandler();
    preco->SetExcitationHandler(handler);
    
    // Nuclear evaporation
    handler->SetEvaporation(new G4Evaporation());
    handler->SetFermiModel(new G4FermiBreakUpModel());
    handler->SetMultiFragmentation(new G4StatMF());
}
```

### Atomic Scale Physics APIs (MeV → keV)

#### Electromagnetic Process Implementation

```cpp
#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4GammaConversion.hh"
#include "G4RayleighScattering.hh"

// Low-energy EM physics with detailed atomic relaxation
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4PenelopePhotoElectricModel.hh"

void ConstructEMProcesses() {
    auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    
    while((*particleIterator)()) {
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        
        if(particle == G4Gamma::Gamma()) {
            // Photoelectric effect with atomic deexcitation
            G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
            pe->SetEmModel(new G4LivermorePhotoElectricModel());
            pmanager->AddDiscreteProcess(pe);
            
            // Compton scattering
            G4ComptonScattering* cs = new G4ComptonScattering();
            cs->SetEmModel(new G4LivermoreComptonModel());
            pmanager->AddDiscreteProcess(cs);
            
            // Pair production
            G4GammaConversion* gc = new G4GammaConversion();
            pmanager->AddDiscreteProcess(gc);
            
            // Rayleigh scattering
            G4RayleighScattering* rs = new G4RayleighScattering();
            rs->SetEmModel(new G4LivermoreRayleighModel());
            pmanager->AddDiscreteProcess(rs);
        }
    }
}
```

#### Atomic Deexcitation API

```cpp
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

// Enable atomic deexcitation globally
void EnableAtomicDeexcitation() {
    G4UAtomicDeexcitation* deex = new G4UAtomicDeexcitation();
    deex->SetFluo(true);        // Fluorescence
    deex->SetAuger(true);       // Auger electrons
    deex->SetPIXE(true);        // Particle induced X-ray emission
    
    G4LossTableManager::Instance()->SetAtomDeexcitation(deex);
}

// Region-specific activation
void SetDeexcitationForRegion(G4Region* region) {
    G4ProductionCutsTable::GetProductionCutsTable()
        ->GetDefaultRegion()
        ->GetDeexcitationProcesses()
        ->SetFluo(true);
}
```

## Atomic Particles to Molecular Systems - APIs

### Material Definition APIs

```cpp
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

// Method 1: Using NIST database
G4NistManager* nist = G4NistManager::Instance();
G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
G4Material* dna = nist->FindOrBuildMaterial("G4_DNA");

// Method 2: Custom molecular material
G4Material* CreateCustomMolecule() {
    G4double z, a;
    G4String name, symbol;
    
    // Define elements
    G4Element* H = new G4Element(name="Hydrogen", symbol="H", z=1., a=1.008*g/mole);
    G4Element* C = new G4Element(name="Carbon", symbol="C", z=6., a=12.01*g/mole);
    G4Element* O = new G4Element(name="Oxygen", symbol="O", z=8., a=16.00*g/mole);
    G4Element* N = new G4Element(name="Nitrogen", symbol="N", z=7., a=14.01*g/mole);
    
    // Create glucose (C6H12O6)
    G4Material* glucose = new G4Material(name="Glucose", 1.54*g/cm3, 3);
    glucose->AddElement(C, 6);
    glucose->AddElement(H, 12);
    glucose->AddElement(O, 6);
    
    return glucose;
}

// Method 3: Material with custom mean excitation energy
void SetMolecularProperties(G4Material* material) {
    G4double meanExcitationEnergy = 75.0*eV;  // Custom I-value
    material->GetIonisation()->SetMeanExcitationEnergy(meanExcitationEnergy);
    
    // Set density effect parameters
    G4DensityEffectData* densityEffectData = 
        G4NistManager::Instance()->GetDensityEffectData();
    densityEffectData->AddMaterial(material);
}
```

### Geant4-DNA Extension APIs

```cpp
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAElastic.hh"
#include "G4DNAExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"

// Initialize DNA physics
class G4EmDNAPhysics : public G4VPhysicsConstructor {
public:
    void ConstructProcess() override {
        G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
        auto particleIterator = GetParticleIterator();
        particleIterator->reset();
        
        while((*particleIterator)()) {
            G4ParticleDefinition* particle = particleIterator->value();
            G4String particleName = particle->GetParticleName();
            
            if(particleName == "e-") {
                // DNA-specific electron processes
                ph->RegisterProcess(new G4DNAElastic("e-_G4DNAElastic"), particle);
                ph->RegisterProcess(new G4DNAExcitation("e-_G4DNAExcitation"), particle);
                ph->RegisterProcess(new G4DNAIonisation("e-_G4DNAIonisation"), particle);
                ph->RegisterProcess(new G4DNAElectronSolvation("e-_G4DNASolvation"), particle);
                ph->RegisterProcess(new G4DNAVibExcitation("e-_G4DNAVibExcitation"), particle);
            }
            else if(particleName == "proton") {
                ph->RegisterProcess(new G4DNAExcitation("proton_G4DNAExcitation"), particle);
                ph->RegisterProcess(new G4DNAIonisation("proton_G4DNAIonisation"), particle);
                ph->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), particle);
            }
        }
    }
};

// Chemistry simulation
void InitializeDNAChemistry() {
    G4DNAChemistryManager* chemMan = G4DNAChemistryManager::Instance();
    chemMan->SetChemistryList(new G4EmDNAChemistry());
    chemMan->Initialize();
    
    // Define chemical reactions
    G4DNAMolecularReactionTable* reactionTable = 
        G4DNAMolecularReactionTable::GetReactionTable();
    
    // OH• + OH• → H2O2
    reactionTable->SetReaction(
        G4MoleculeTable::Instance()->GetConfiguration("OH"),
        G4MoleculeTable::Instance()->GetConfiguration("OH"),
        2.5e10 * (1e-3*m3/(mole*s)),  // Reaction rate
        G4MoleculeTable::Instance()->GetConfiguration("H2O2")
    );
}
```

### Optical Physics for Molecular Spectroscopy

```cpp
#include "G4OpticalPhysics.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpWLS.hh"
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"

// Configure optical properties
class OpticalPhysics : public G4OpticalPhysics {
public:
    OpticalPhysics() : G4OpticalPhysics() {
        SetMaxNumPhotonsPerStep(100);
        SetMaxBetaChangePerStep(10.0);
        SetTrackSecondariesFirst(kCerenkov, true);
        SetTrackSecondariesFirst(kScintillation, true);
        
        // Configure wavelength shifting (fluorescence)
        SetWLSTimeProfile("exponential");
        SetScintillationYieldFactor(1.0);
        SetScintillationExcitationRatio(0.0);
        
        // Enable specific processes
        Configure(kWLS, true);
        Configure(kCerenkov, true);
        Configure(kScintillation, true);
        Configure(kAbsorption, true);
        Configure(kRayleigh, true);
        Configure(kMieHG, false);
        Configure(kBoundary, true);
    }
};

// Material optical properties
void SetOpticalProperties(G4Material* material) {
    const G4int nEntries = 50;
    G4double photonEnergy[nEntries];
    G4double refractiveIndex[nEntries];
    G4double absorption[nEntries];
    G4double scintillation[nEntries];
    
    for(G4int i = 0; i < nEntries; i++) {
        photonEnergy[i] = (1.0 + i * 0.1) * eV;
        refractiveIndex[i] = 1.58;
        absorption[i] = 50.0 * cm;
        scintillation[i] = 1.0;
    }
    
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries);
    mpt->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries);
    mpt->AddProperty("FASTCOMPONENT", photonEnergy, scintillation, nEntries);
    mpt->AddConstProperty("SCINTILLATIONYIELD", 10000.0/MeV);
    mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    mpt->AddConstProperty("FASTTIMECONSTANT", 1.0*ns);
    
    material->SetMaterialPropertiesTable(mpt);
}
```

## Interaction Physics Implementation

### Electromagnetic Interaction APIs

#### Photon Process Registration

```cpp
#include "G4EmBuilder.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

// Choose EM physics constructor based on energy range
void SelectEMPhysics(G4VModularPhysicsList* physicsList, const G4String& option) {
    if(option == "standard") {
        physicsList->RegisterPhysics(new G4EmStandardPhysics());
    }
    else if(option == "livermore") {
        // Best for low energies (250 eV - 100 GeV)
        physicsList->RegisterPhysics(new G4EmLivermorePhysics());
    }
    else if(option == "penelope") {
        // Detailed atomic relaxation (100 eV - 1 GeV)
        physicsList->RegisterPhysics(new G4EmPenelopePhysics());
    }
    else if(option == "option4") {
        // Most accurate, slowest
        physicsList->RegisterPhysics(new G4EmStandardPhysics_option4());
    }
}

// Custom EM process configuration
#include "G4EmParameters.hh"

void ConfigureEMParameters() {
    G4EmParameters* emParams = G4EmParameters::Instance();
    
    // Global settings
    emParams->SetMinEnergy(100*eV);
    emParams->SetMaxEnergy(10*TeV);
    emParams->SetNumberOfBinsPerDecade(20);
    
    // Multiple scattering
    emParams->SetMscRangeFactor(0.04);
    emParams->SetMscStepLimitType(fUseDistanceToBoundary);
    
    // Ionization and atomic deexcitation
    emParams->SetFluo(true);
    emParams->SetAuger(true);
    emParams->SetPIXE(true);
    emParams->SetDeexcitationIgnoreCut(true);
    
    // Energy loss
    emParams->SetLossFluctuations(true);
    emParams->SetUseCutAsFinalRange(false);
    emParams->SetLinearLossLimit(0.01);
}
```

### Hadronic Interaction APIs

#### Neutron Physics Implementation

```cpp
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPThermalScattering.hh"

// High-precision neutron physics (thermal to 20 MeV)
class NeutronHPPhysics : public G4VPhysicsConstructor {
public:
    void ConstructProcess() override {
        G4ProcessManager* pManager = G4Neutron::Neutron()->GetProcessManager();
        
        // Elastic scattering with S(α,β) for molecular systems
        G4NeutronHPElastic* elasticHP = new G4NeutronHPElastic();
        G4NeutronHPElasticData* elasticData = new G4NeutronHPElasticData();
        elasticHP->AddDataSet(elasticData);
        
        // Thermal scattering for molecular moderators
        G4NeutronHPThermalScattering* thermal = new G4NeutronHPThermalScattering();
        G4NeutronHPThermalScatteringData* thermalData = 
            new G4NeutronHPThermalScatteringData();
        thermal->AddDataSet(thermalData);
        
        // Inelastic
        G4NeutronHPInelastic* inelasticHP = new G4NeutronHPInelastic();
        
        // Capture
        G4NeutronHPCapture* captureHP = new G4NeutronHPCapture();
        
        // Fission (if applicable)
        G4NeutronHPFission* fissionHP = new G4NeutronHPFission();
        
        // Register processes
        G4HadronElasticProcess* elasticProcess = new G4HadronElasticProcess();
        elasticProcess->RegisterMe(elasticHP);
        pManager->AddDiscreteProcess(elasticProcess);
        
        G4NeutronInelasticProcess* inelasticProcess = new G4NeutronInelasticProcess();
        inelasticProcess->RegisterMe(inelasticHP);
        pManager->AddDiscreteProcess(inelasticProcess);
        
        G4HadronCaptureProcess* captureProcess = new G4HadronCaptureProcess();
        captureProcess->RegisterMe(captureHP);
        pManager->AddDiscreteProcess(captureProcess);
    }
};

// Data library configuration
void SetupNeutronData() {
    // Point to G4NDL data
    setenv("G4NEUTRONHPDATA", "/path/to/G4NDL4.6", 1);
    
    // Use specific evaluation
    G4NeutronHPManager::GetInstance()->SetUseWendtFissionModel(true);
    G4NeutronHPManager::GetInstance()->SetProduceFissionFragments(true);
    G4NeutronHPManager::GetInstance()->SetSkipMissingIsotopes(false);
    G4NeutronHPManager::GetInstance()->SetDoNotAdjustFinalState(false);
    G4NeutronHPManager::GetInstance()->SetUseOnlyPhotoEvaporation(false);
    G4NeutronHPManager::GetInstance()->SetNeglectDoppler(false);
}
```

#### Ion and Proton Interactions

```cpp
#include "G4IonPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonQMDPhysics.hh"

// Ion physics with QMD for molecular targets
class IonPhysics : public G4VPhysicsConstructor {
public:
    void ConstructProcess() override {
        G4ProcessManager* pManager = nullptr;
        
        // Add processes for all ions
        auto particleIterator = GetParticleIterator();
        particleIterator->reset();
        
        while((*particleIterator)()) {
            G4ParticleDefinition* particle = particleIterator->value();
            
            if(particle->GetAtomicMass() > 1) {  // Ions heavier than proton
                pManager = particle->GetProcessManager();
                
                // Elastic
                G4HadronElasticProcess* elasticProcess = new G4HadronElasticProcess();
                elasticProcess->RegisterMe(new G4IonElasticModel());
                pManager->AddDiscreteProcess(elasticProcess);
                
                // Inelastic - QMD for light ions
                G4IonInelasticProcess* inelasticProcess = new G4IonInelasticProcess();
                
                if(particle->GetAtomicMass() < 18) {
                    // Use QMD for light ions (better for molecular targets)
                    G4QMDReaction* qmd = new G4QMDReaction();
                    qmd->SetMinEnergy(0.0*MeV);
                    qmd->SetMaxEnergy(10.0*GeV);
                    inelasticProcess->RegisterMe(qmd);
                } else {
                    // Use Binary cascade for heavy ions
                    G4BinaryLightIonReaction* binary = new G4BinaryLightIonReaction();
                    binary->SetMinEnergy(0.0*MeV);
                    binary->SetMaxEnergy(10.0*GeV);
                    inelasticProcess->RegisterMe(binary);
                }
                
                pManager->AddDiscreteProcess(inelasticProcess);
            }
        }
    }
};
```

## Physics Lists and Model APIs

### Modular Physics List Construction

```cpp
#include "G4VModularPhysicsList.hh"
#include "G4PhysListFactory.hh"
#include "G4PhysicsConstructorRegistry.hh"

// Factory method for standard physics lists
G4VModularPhysicsList* GetPhysicsList(const G4String& name) {
    G4PhysListFactory factory;
    
    // Available: FTFP_BERT, QGSP_BERT, QGSP_BIC_HP, Shielding, etc.
    if(factory.IsReferencePhysList(name)) {
        return factory.GetReferencePhysList(name);
    }
    
    // Custom physics list
    return new MyCustomPhysicsList();
}

// Building custom modular physics list
class MyCustomPhysicsList : public G4VModularPhysicsList {
public:
    MyCustomPhysicsList() {
        defaultCutValue = 0.7*mm;
        
        // Register physics constructors
        RegisterPhysics(new G4EmStandardPhysics_option4());
        RegisterPhysics(new G4EmExtraPhysics());          // Synchrotron, muon-nuclear
        RegisterPhysics(new G4OpticalPhysics());
        RegisterPhysics(new G4HadronElasticPhysicsHP());
        RegisterPhysics(new G4HadronPhysicsFTFP_BERT_HP());
        RegisterPhysics(new G4IonBinaryCascadePhysics());
        RegisterPhysics(new G4NeutronTrackingCut());
        RegisterPhysics(new G4DecayPhysics());
        RegisterPhysics(new G4RadioactiveDecayPhysics());
        RegisterPhysics(new G4StoppingPhysics());         // Muon/hadron stopping
        
        // DNA physics for molecular damage
        RegisterPhysics(new G4EmDNAPhysics());
    }
    
    virtual void SetCuts() override {
        SetCutValue(defaultCutValue, "gamma");
        SetCutValue(defaultCutValue, "e-");
        SetCutValue(defaultCutValue, "e+");
        SetCutValue(defaultCutValue, "proton");
        
        // Region-specific cuts
        G4Region* region = G4RegionStore::GetInstance()->GetRegion("TargetRegion");
        if(region) {
            G4ProductionCuts* cuts = new G4ProductionCuts();
            cuts->SetProductionCut(0.01*mm, G4ProductionCuts::GetIndex("gamma"));
            cuts->SetProductionCut(0.01*mm, G4ProductionCuts::GetIndex("e-"));
            region->SetProductionCuts(cuts);
        }
    }
};
```

### Physics Process Management

```cpp
#include "G4ProcessTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

// Runtime physics configuration
void ModifyPhysicsProcesses() {
    G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();
    
    // Find and modify specific process
    G4VProcess* process = processTable->FindProcess("phot", G4Gamma::Gamma());
    if(process) {
        G4PhotoElectricEffect* photoEffect = 
            dynamic_cast<G4PhotoElectricEffect*>(process);
        if(photoEffect) {
            // Switch to Penelope model for better low-energy accuracy
            photoEffect->SetEmModel(new G4PenelopePhotoElectricModel());
        }
    }
    
    // Deactivate process for specific particle
    G4ProcessManager* pManager = G4Electron::Electron()->GetProcessManager();
    G4ProcessVector* processVector = pManager->GetProcessList();
    for(size_t i = 0; i < processVector->size(); i++) {
        G4VProcess* proc = (*processVector)[i];
        if(proc->GetProcessName() == "eBrem") {
            pManager->SetProcessActivation(proc, false);
        }
    }
}

// Process verbosity control
void SetProcessVerbosity() {
    G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();
    processTable->SetVerboseLevel(1);
    
    // Individual process verbosity
    G4VProcess* compton = processTable->FindProcess("compt", G4Gamma::Gamma());
    if(compton) compton->SetVerboseLevel(2);
}
```

## Bindings and Extensions

### Python Bindings (Geant4Py)

```python
# Python interface to Geant4
import Geant4 as g4
from Geant4 import mm, cm, m, MeV, GeV

# Create run manager
runManager = g4.G4RunManager()

# Detector construction in Python
class MyDetector(g4.G4VUserDetectorConstruction):
    def __init__(self):
        g4.G4VUserDetectorConstruction.__init__(self)
        
    def Construct(self):
        # Materials
        nist = g4.G4NistManager.Instance()
        water = nist.FindOrBuildMaterial("G4_WATER")
        
        # Geometry
        worldBox = g4.G4Box("World", 1*m, 1*m, 1*m)
        worldLV = g4.G4LogicalVolume(worldBox, water, "World")
        worldPV = g4.G4PVPlacement(None, g4.G4ThreeVector(), 
                                  worldLV, "World", None, False, 0)
        return worldPV

# Physics list
physicsList = g4.FTFP_BERT()
runManager.SetUserInitialization(physicsList)

# Primary generator
class MyPrimaryGenerator(g4.G4VUserPrimaryGeneratorAction):
    def __init__(self):
        g4.G4VUserPrimaryGeneratorAction.__init__(self)
        self.particleGun = g4.G4ParticleGun(1)
        
    def GeneratePrimaries(self, event):
        self.particleGun.SetParticleDefinition(g4.G4Gamma.Gamma())
        self.particleGun.SetParticleEnergy(1*MeV)
        self.particleGun.SetParticlePosition(g4.G4ThreeVector(0,0,0))
        self.particleGun.GeneratePrimaryVertex(event)

# Run simulation
runManager.SetUserInitialization(MyDetector())
runManager.SetUserAction(MyPrimaryGenerator())
runManager.Initialize()
runManager.BeamOn(1000)
```

### ROOT Integration

```cpp
#include "g4root.hh"  // or g4xml.hh, g4csv.hh, g4hdf5.hh

// Analysis manager for data output
class RunAction : public G4UserRunAction {
private:
    G4AnalysisManager* analysisManager;
    
public:
    RunAction() {
        analysisManager = G4AnalysisManager::Instance();
        analysisManager->SetVerboseLevel(1);
        analysisManager->SetNtupleMerging(true);
        
        // Create histograms
        analysisManager->CreateH1("EdepHist", "Energy deposition", 100, 0., 100*MeV);
        analysisManager->CreateH2("XYDist", "XY distribution", 
                                  100, -10*cm, 10*cm, 100, -10*cm, 10*cm);
        
        // Create ntuple
        analysisManager->CreateNtuple("SimData", "Simulation data");
        analysisManager->CreateNtupleDColumn("Energy");
        analysisManager->CreateNtupleDColumn("X");
        analysisManager->CreateNtupleDColumn("Y");
        analysisManager->CreateNtupleDColumn("Z");
        analysisManager->FinishNtuple();
    }
    
    void BeginOfRunAction(const G4Run*) override {
        analysisManager->OpenFile("output.root");
    }
    
    void EndOfRunAction(const G4Run*) override {
        analysisManager->Write();
        analysisManager->CloseFile();
    }
};

// Fill histograms in stepping action
void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4double edep = step->GetTotalEnergyDeposit();
    if(edep > 0) {
        auto analysisManager = G4AnalysisManager::Instance();
        analysisManager->FillH1(0, edep);
        
        G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
        analysisManager->FillH2(0, pos.x(), pos.y());
        
        analysisManager->FillNtupleDColumn(0, edep);
        analysisManager->FillNtupleDColumn(1, pos.x());
        analysisManager->FillNtupleDColumn(2, pos.y());
        analysisManager->FillNtupleDColumn(3, pos.z());
        analysisManager->AddNtupleRow();
    }
}
```

### GDML (Geometry Description Markup Language)

```cpp
#include "G4GDMLParser.hh"

// Export geometry to GDML
void ExportGeometry(G4VPhysicalVolume* world) {
    G4GDMLParser parser;
    parser.Write("geometry.gdml", world);
}

// Import geometry from GDML
G4VPhysicalVolume* ImportGeometry() {
    G4GDMLParser parser;
    parser.Read("geometry.gdml");
    return parser.GetWorldVolume();
}
```

```xml
<!-- Example GDML file -->
<gdml>
  <define>
    <constant name="world_size" value="1000"/>
  </define>
  
  <materials>
    <material name="Water" Z="1">
      <D value="1.0"/>
      <atom value="18.0"/>
    </material>
  </materials>
  
  <solids>
    <box name="World" x="world_size" y="world_size" z="world_size"/>
    <sphere name="Target" rmax="100"/>
  </solids>
  
  <structure>
    <volume name="Target_log">
      <materialref ref="Water"/>
      <solidref ref="Target"/>
    </volume>
    
    <volume name="World_log">
      <materialref ref="Air"/>
      <solidref ref="World"/>
      <physvol>
        <volumeref ref="Target_log"/>
        <position x="0" y="0" z="0"/>
      </physvol>
    </volume>
  </structure>
  
  <setup name="Default" version="1.0">
    <world ref="World_log"/>
  </setup>
</gdml>
```

### Visualization Extensions

```cpp
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"

// Visualization setup
void SetupVisualization(int argc, char** argv) {
    // Visualization manager
    G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
    
    // UI manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    // OpenGL visualization
    UImanager->ApplyCommand("/vis/open OGL 800x600-0+0");
    UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 90 180");
    UImanager->ApplyCommand("/vis/drawVolume");
    UImanager->ApplyCommand("/vis/viewer/set/style wireframe");
    UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
    UImanager->ApplyCommand("/vis/scene/add/hits");
    
    // HepRep for offline viewing
    UImanager->ApplyCommand("/vis/open HepRepFile");
    UImanager->ApplyCommand("/vis/heprep/setCullInvisibles true");
    
    // VRML2 export
    UImanager->ApplyCommand("/vis/open VRML2FILE");
    UImanager->ApplyCommand("/vis/viewer/flush");
}

// Qt-based UI
void RunWithQt(int argc, char** argv) {
    G4UIExecutive* ui = new G4UIExecutive(argc, argv, "Qt");
    
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/control/execute init.mac");
    
    ui->SessionStart();
    delete ui;
}
```

## Advanced API Usage and Custom Physics

### Custom Physics Process

```cpp
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"

// Custom molecular excitation process
class MolecularExcitation : public G4VDiscreteProcess {
public:
    MolecularExcitation(const G4String& name = "MolecularExcitation")
        : G4VDiscreteProcess(name) {}
    
    G4double GetMeanFreePath(const G4Track& track, 
                            G4double, 
                            G4ForceCondition* condition) override {
        // Calculate mean free path based on cross-section
        G4double energy = track.GetKineticEnergy();
        G4double crossSection = CalculateCrossSection(energy);
        G4double density = track.GetMaterial()->GetDensity();
        
        return crossSection > 0 ? 1./(density * crossSection) : DBL_MAX;
    }
    
    G4VParticleChange* PostStepDoIt(const G4Track& track, 
                                   const G4Step& step) override {
        aParticleChange.Initialize(track);
        
        // Implement excitation physics
        G4double excitationEnergy = 5.0*eV;  // Example
        G4double newKineticEnergy = track.GetKineticEnergy() - excitationEnergy;
        
        if(newKineticEnergy > 0) {
            aParticleChange.ProposeEnergy(newKineticEnergy);
            
            // Create de-excitation products (e.g., photon)
            G4DynamicParticle* photon = new G4DynamicParticle(
                G4Gamma::Gamma(),
                G4RandomDirection(),
                excitationEnergy
            );
            aParticleChange.AddSecondary(photon);
        } else {
            aParticleChange.ProposeTrackStatus(fStopAndKill);
        }
        
        return &aParticleChange;
    }
    
private:
    G4double CalculateCrossSection(G4double energy) {
        // Implement cross-section calculation
        return 1.0e-16 * cm2 * exp(-energy/10.0*eV);
    }
};

// Register custom process
void RegisterCustomProcess() {
    G4ProcessManager* pManager = G4Electron::Electron()->GetProcessManager();
    pManager->AddDiscreteProcess(new MolecularExcitation());
}
```

### Biasing and Variance Reduction

```cpp
#include "G4VBiasingOperator.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4GenericBiasingPhysics.hh"

// Importance sampling for molecular regions
class MolecularBiasingOperator : public G4VBiasingOperator {
public:
    MolecularBiasingOperator() : G4VBiasingOperator("MolecularBiasing") {}
    
    G4VBiasingOperation* ProposeOccurenceBiasingOperation(
        const G4Track* track,
        const G4BiasingProcessInterface* callingProcess) override {
        
        // Enhance specific processes in molecular region
        if(track->GetVolume()->GetName() == "MolecularTarget") {
            G4String processName = callingProcess->GetWrappedProcess()->GetProcessName();
            
            if(processName == "compt" || processName == "phot") {
                // Increase interaction probability by factor of 10
                return new G4BOptnForceCollision(callingProcess->GetWrappedProcess(), 10.0);
            }
        }
        return nullptr;
    }
};

// Apply biasing to physics list
void AddBiasing(G4VModularPhysicsList* physicsList) {
    G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
    biasingPhysics->BiasAllNeutral();
    biasingPhysics->BiasAllCharged();
    physicsList->RegisterPhysics(biasingPhysics);
}
```

### Parallel Computing APIs

```cpp
#include "G4MTRunManager.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

// Thread-safe singleton for shared data
class SharedData {
private:
    static SharedData* instance;
    static G4Mutex instanceMutex;
    std::vector<G4double> data;
    G4Mutex dataMutex;
    
public:
    static SharedData* GetInstance() {
        G4AutoLock lock(&instanceMutex);
        if(!instance) instance = new SharedData();
        return instance;
    }
    
    void AddData(G4double value) {
        G4AutoLock lock(&dataMutex);
        data.push_back(value);
    }
};

// Multi-threaded run manager setup
void SetupMT() {
    G4MTRunManager* runManager = new G4MTRunManager;
    
    // Set number of threads
    G4int nThreads = G4Threading::G4GetNumberOfCores();
    runManager->SetNumberOfThreads(nThreads);
    
    // Thread-specific random seeds
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/random/setSeeds 12345 67890");
}
```

### GPU Acceleration (CUDA Extension)

```cpp
// G4CudaInterface.hh - Custom CUDA integration
#ifdef G4CUDA_ENABLED
#include <cuda_runtime.h>

class G4CudaManager {
public:
    static void InitializeCuda() {
        int deviceCount;
        cudaGetDeviceCount(&deviceCount);
        
        if(deviceCount > 0) {
            cudaSetDevice(0);
            G4cout << "CUDA initialized with " << deviceCount << " devices" << G4endl;
        }
    }
    
    static void LaunchKernel(/* kernel parameters */) {
        // Launch CUDA kernel for parallel processing
        // e.g., cross-section calculations, ray tracing
    }
};
#endif
```

### Machine Learning Integration

```cpp
// Integration with ML frameworks for physics modeling
#include "G4MLPhysicsInterface.hh"  // Hypothetical ML interface

class MLCrossSection : public G4VCrossSectionDataSet {
private:
    void* model;  // TensorFlow or PyTorch model
    
public:
    MLCrossSection() {
        // Load trained model
        model = LoadModel("molecular_xsection.pb");
    }
    
    G4double GetElementCrossSection(const G4DynamicParticle* particle,
                                   G4int Z,
                                   const G4Material* mat) override {
        // Use ML model to predict cross-section
        G4double energy = particle->GetKineticEnergy();
        G4double prediction = EvaluateModel(model, energy, Z, mat->GetDensity());
        return prediction * barn;
    }
};
```

### Web-Based Interfaces

```cpp
// REST API for remote simulation control
#include "G4RestServer.hh"  // Custom REST server

class SimulationServer {
public:
    void StartServer(int port = 8080) {
        // Setup HTTP endpoints
        server.Get("/status", [](const Request& req, Response& res) {
            res.set_content(GetSimulationStatus(), "application/json");
        });
        
        server.Post("/run", [](const Request& req, Response& res) {
            G4UImanager::GetUIpointer()->ApplyCommand("/run/beamOn 1000");
            res.set_content("Simulation started", "text/plain");
        });
        
        server.listen("localhost", port);
    }
};
```

## Extension Frameworks

### GATE (Medical Physics)
```cpp
// GATE extensions for medical imaging
#include "GateConfiguration.h"
#include "GateSourceMgr.hh"
#include "GateDetectorConstruction.hh"

// PET scanner simulation
GateConfiguration* gate = new GateConfiguration();
gate->SetupPETScanner();
gate->SetActivity(10*MBq);
gate->SetAcquisitionTime(600*s);
```

### TOPAS (Particle Therapy)
```cpp
// TOPAS parameter system
// Parameters defined in text files:
// s:Ge/World/Material = "Air"
// d:Ge/World/HLX = 2.0 m
// i:Ts/NumberOfThreads = 8
```

### Neutrino and Astroparticle Extensions
```cpp
// GLG4sim for neutrino detectors
#include "GLG4Scint.hh"
#include "GLG4OpticalPhysics.hh"

// WCSim for water Cherenkov detectors
#include "WCSimDetectorConstruction.hh"
#include "WCSimPrimaryGeneratorAction.hh"
```

---

This comprehensive guide demonstrates how Geant4's extensive API system enables simulation across all scales, from TeV accelerator physics to eV molecular interactions, with bindings and extensions for specialized applications.

I've extensively expanded the markdown guide to include direct correlations with Geant4 APIs, showing exactly how each physics concept maps to actual C++ implementation. The document now covers:

## Key API Additions:

### 1. **Core Architecture Implementation**
- Specific manager classes (`G4RunManager`, `G4MTRunManager`, etc.)
- Mandatory user action classes with code examples
- Proper initialization patterns

### 2. **Physics Scale Implementations**
- **High Energy**: String model APIs (`G4FTFModel`, `G4QGSMFragmentation`)
- **Cascade Models**: Bertini and Binary cascade setup
- **Atomic Scale**: Livermore/Penelope model configuration
- **Molecular**: DNA extension APIs and chemistry managers

### 3. **Material and Molecular System APIs**
- NIST material database usage
- Custom molecular material creation
- Optical properties configuration
- DNA-specific process registration

### 4. **Interaction Physics with Code**
- Electromagnetic process registration
- Neutron HP (High Precision) setup
- Ion and QMD physics implementation
- Complete physics list construction

### 5. **Bindings and Extensions**
- **Python**: Geant4Py examples with complete simulation
- **ROOT**: Analysis manager integration
- **GDML**: Geometry import/export
- **Visualization**: OpenGL, Qt, HepRep setup

### 6. **Advanced Features**
- Custom physics process development
- Biasing and variance reduction
- Multi-threading configuration
- GPU acceleration hooks
- Machine learning integration concepts
- Web API interfaces

### 7. **Specialized Frameworks**
- GATE for medical physics
- TOPAS for particle therapy
- GLG4sim/WCSim for neutrino physics

Each section now includes:
- Relevant `#include` statements
- Class inheritance patterns
- Method implementations
- Configuration examples
- Integration patterns

This provides developers with a complete reference showing not just what Geant4 can simulate, but exactly how to implement it in code, from the highest energy particle collisions down to molecular-level interactions.