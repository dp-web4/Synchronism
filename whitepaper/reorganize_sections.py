#!/usr/bin/env python3
"""
Reorganize Synchronism sections into fractal directory structure
"""

import os
import shutil
import re

# Define the fractal structure
STRUCTURE = {
    "00-executive-summary": {
        "files": ["00-executive-summary.md"],
        "subdirs": {}
    },
    "01-introduction": {
        "files": ["01-introduction.md", "02-introduction.md"],
        "subdirs": {}
    },
    "02-foundations": {
        "files": [],
        "subdirs": {
            "01-perspective": {
                "files": ["03-contents.md"],
                "subdirs": {}
            },
            "02-hermetic-principles": {
                "files": ["03-hermetic-principles.md"],
                "subdirs": {}
            }
        }
    },
    "03-core-concepts": {
        "files": [],
        "subdirs": {
            "01-intent": {
                "files": ["04-intent-quantization.md", "05-intent-saturation.md"],
                "subdirs": {}
            },
            "02-coherence": {
                "files": ["07-7-coherence-and-feedback.md", "89-coherence-function.md"],
                "subdirs": {}
            },
            "03-markov-systems": {
                "files": [
                    "08-8-markov-blankets-and-scale-bo.md",
                    "09-example-markov-blankets-in-a-n.md",
                    "10-9-markov-relevancy-horizon.md",
                    "11-example-markov-relevancy-horiz.md"
                ],
                "subdirs": {}
            },
            "04-spectral-existence": {
                "files": [
                    "12-10-spectral-existence-in-synch.md",
                    "13-existence-as-a-spectrum.md",
                    "14-dark-matter-and-dark-energy-as.md",
                    "15-persistence-and-decoherence.md",
                    "19-spectral-existence-and-perceiv.md"
                ],
                "subdirs": {}
            },
            "05-abstraction": {
                "files": ["16-11-abstraction.md"],
                "subdirs": {}
            }
        }
    },
    "04-quantum-perspectives": {
        "files": ["07-quantum-cosmic-bridge.md"],
        "subdirs": {
            "01-quantum-phenomena": {
                "files": [
                    "20-superposition-as-spectral-exis.md",
                    "21-wave-particle-duality-and-spec.md",
                    "22-entanglement-as-synchronizatio.md",
                    "23-crt-entanglement-analogy.md"
                ],
                "subdirs": {}
            },
            "02-witness-effect": {
                "files": [
                    "24-witness-synchronization-and-pe.md",
                    "25-5-witness-effect.md"
                ],
                "subdirs": {}
            },
            "03-relativity": {
                "files": [
                    "26-speed-limits-and-time-dilation.md",
                    "27-pendulum-clock-analogy.md",
                    "28-speed-of-light-and-time-dilati.md",
                    "29-velocity-and-complexity.md"
                ],
                "subdirs": {}
            }
        }
    },
    "05-macro-phenomena": {
        "files": ["17-and-macro-phenomena.md"],
        "subdirs": {
            "01-decoherence": {
                "files": [
                    "35-macro-decoherence-across-scale.md",
                    "36-macro-decoherence-and-high-spe.md"
                ],
                "subdirs": {}
            },
            "02-phase-transitions": {
                "files": [
                    "37-temperature-as-speed-of-intent.md",
                    "38-phase-transitions-as-macro-coh.md",
                    "39-emergent-group-behavior-in-pha.md",
                    "75-5192-phase-transitions.md",
                    "76-5193-critical-phenomena.md"
                ],
                "subdirs": {}
            },
            "03-energy": {
                "files": [
                    "40-10-energy-in-synchronism.md",
                    "41-energy-as-magnitude-of-intent-.md",
                    "42-localized-and-resonant-energy.md",
                    "43-expanding-on-energy-types.md",
                    "44-bulk-energy-and-mass.md",
                    "45-conversion-between-localized-a.md",
                    "68-refinement-on-energy.md",
                    "69-5181-fundamental-definition.md",
                    "70-5182-forms-of-energy.md",
                    "71-5183-energy-conservation-and-t.md",
                    "72-5184-entropy-emergence-and-the.md"
                ],
                "subdirs": {}
            }
        }
    },
    "06-applications": {
        "files": ["07-implications.md"],
        "subdirs": {
            "01-chemistry": {
                "files": [
                    "47-2-chemistry-in-synchronism.md",
                    "48-chemical-bonding-and-molecular.md",
                    "49-reinterpreting-chemical-reacti.md",
                    "50-chemical-reactions.md",
                    "51-catalysis.md",
                    "52-markov-relevancy-horizon-in-ch.md"
                ],
                "subdirs": {}
            },
            "02-biology": {
                "files": [
                    "53-13-coherence-of-life-and-cogni.md",
                    "54-temperature-energy-and-biologi.md",
                    "55-cognitive-coherence-and-its-fr.md",
                    "77-refinement-on-coherence-of-lif.md",
                    "78-5201-biological-coherence.md",
                    "79-5202-cognitive-coherence.md",
                    "80-5203-resilience-and-adaptation.md"
                ],
                "subdirs": {}
            },
            "03-cosmology": {
                "files": [
                    "33-cosmological-implications.md",
                    "56-impact-of-high-speed-travel-on.md",
                    "57-black-holes-and-dark-matter-gr.md",
                    "59-gravitational-lensing-anomalie.md",
                    "60-galaxy-rotation-curves.md",
                    "61-microlensing-events.md",
                    "62-supermassive-black-holes-and-d.md",
                    "63-cosmic-microwave-background-cm.md"
                ],
                "subdirs": {}
            },
            "04-technology": {
                "files": [
                    "30-572-applications-and-implicati.md",
                    "31-high-speed-travel-and-space-ex.md",
                    "32-advanced-computational-models.md",
                    "64-superconductivity-in-synchroni.md"
                ],
                "subdirs": {}
            },
            "05-electromagnetism": {
                "files": [
                    "65-5171-core-principles.md",
                    "66-5172-reinterpreting-maxwells-e.md",
                    "67-5173-electromagnetic-interacti.md",
                    "46-universal-field-in-synchronism.md"
                ],
                "subdirs": {}
            }
        }
    },
    "07-mathematical-framework": {
        "files": ["08-mathematical-framework.md"],
        "subdirs": {
            "01-core-formalism": {
                "files": [
                    "82-proposed-mathematical-formalis.md",
                    "83-tensor-representation-of-tensi.md",
                    "84-local-transfer-potential.md",
                    "85-intent-updating-rule.md"
                ],
                "subdirs": {}
            },
            "02-analysis-methods": {
                "files": [
                    "86-fourier-analysis-for-pattern-s.md",
                    "87-markov-analysis-for-state-tran.md",
                    "88-2-mathematical-representation-.md"
                ],
                "subdirs": {}
            },
            "03-coherence-mathematics": {
                "files": [
                    "90-relationship-to-tension-field.md",
                    "91-example-coherence-in-a-biologi.md",
                    "92-scale-dependent-coherence-func.md",
                    "93-modified-updating-rules.md",
                    "94-coherence-correlation-function.md",
                    "95-order-parameter.md"
                ],
                "subdirs": {}
            },
            "04-spectral-mathematics": {
                "files": [
                    "96-3-mathematical-treatment-of-sp.md",
                    "97-velocity-and-complexity.md",
                    "98-probability-of-transition.md",
                    "99-time-dilation-factor.md",
                    "100-effective-frequency.md",
                    "101-modified-updating-rule.md"
                ],
                "subdirs": {}
            },
            "05-decoherence-mathematics": {
                "files": [
                    "102-4-mathematical-framework-for-m.md",
                    "103-complexity-dependent-decoheren.md",
                    "104-decoherence-probability.md",
                    "105-modification-to-the-coherence-.md",
                    "106-updating-the-intent-field-with.md",
                    "107-effective-time-dilation-with-d.md"
                ],
                "subdirs": {}
            },
            "06-field-theory": {
                "files": [
                    "108-7-tension-field.md",
                    "109-mathematical-treatment-of-grav.md",
                    "110-9-mathematical-treatment-of-su.md",
                    "111-3-components-of-the-interactio.md",
                    "112-4-interaction-tensor-and-spect.md",
                    "113-5-applications-to-dark-matter-.md"
                ],
                "subdirs": {}
            },
            "07-scale-mathematics": {
                "files": [
                    "114-scale-dependent-coherence-matr.md",
                    "115-scale-dependent-feedback-matri.md",
                    "116-emergence-matrix-e.md"
                ],
                "subdirs": {}
            }
        }
    },
    "08-philosophical": {
        "files": ["34-ethical-and-philosophical-cons.md"],
        "subdirs": {}
    },
    "09-conclusion": {
        "files": [
            "11-conclusion.md",
            "120-23-summary-and-integration-of-.md"
        ],
        "subdirs": {}
    },
    "10-future-work": {
        "files": [
            "81-4-open-questions-and-future-di.md",
            "117-20-limitations-and-assumptions.md",
            "118-21-future-directions.md"
        ],
        "subdirs": {}
    }
}

def create_structure(base_path, structure):
    """Recursively create directory structure and move files"""
    for dir_name, content in structure.items():
        dir_path = os.path.join(base_path, dir_name)
        os.makedirs(dir_path, exist_ok=True)
        
        # Move files to this directory
        for filename in content["files"]:
            src = os.path.join("sections", filename)
            if os.path.exists(src):
                dst = os.path.join(dir_path, filename)
                shutil.move(src, dst)
                print(f"Moved {filename} -> {dir_name}/")
        
        # Create subdirectories
        if content["subdirs"]:
            create_structure(dir_path, content["subdirs"])

def create_index_files(base_path):
    """Create index.md files for each directory"""
    for root, dirs, files in os.walk(base_path):
        if root == base_path:
            continue
            
        # Create index for this directory
        index_path = os.path.join(root, "index.md")
        dir_name = os.path.basename(root)
        
        with open(index_path, 'w') as f:
            # Write header
            clean_name = dir_name.split('-', 1)[1] if '-' in dir_name else dir_name
            clean_name = clean_name.replace('-', ' ').title()
            f.write(f"# {clean_name}\n\n")
            
            # List subdirectories
            if dirs:
                f.write("## Sections\n\n")
                for subdir in sorted(dirs):
                    subdir_name = subdir.split('-', 1)[1] if '-' in subdir else subdir
                    subdir_name = subdir_name.replace('-', ' ').title()
                    f.write(f"- [{subdir_name}]({subdir}/index.md)\n")
                f.write("\n")
            
            # List files
            md_files = [f for f in files if f.endswith('.md') and f != 'index.md']
            if md_files:
                f.write("## Documents\n\n")
                for file in sorted(md_files):
                    # Try to extract a better title from the file
                    file_path = os.path.join(root, file)
                    title = file.replace('.md', '').split('-', 1)[1] if '-' in file else file.replace('.md', '')
                    title = title.replace('-', ' ').title()
                    f.write(f"- [{title}]({file})\n")
        
        print(f"Created index: {index_path}")

def main():
    # Create new sections structure
    new_sections = "sections_new"
    os.makedirs(new_sections, exist_ok=True)
    
    # Create the structure
    create_structure(new_sections, STRUCTURE)
    
    # Create index files
    create_index_files(new_sections)
    
    # Move old sections to backup
    if os.path.exists("sections"):
        if os.path.exists("sections_old"):
            shutil.rmtree("sections_old")
        shutil.move("sections", "sections_old")
    
    # Rename new to current
    shutil.move(new_sections, "sections")
    
    print("\n✅ Reorganization complete!")
    print("Old sections backed up to sections_old/")
    print("New fractal structure in sections/")

if __name__ == "__main__":
    main()