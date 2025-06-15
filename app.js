const canvas = document.getElementById('glCanvas');
const container = document.getElementById('simulation-container');

canvas.width = container.clientWidth;
canvas.height = container.clientHeight;

//canvas.width = window.innerWidth;
//canvas.height = window.innerHeight;

const ctx = canvas.getContext('2d');

if (!ctx) {
    throw new Error("Could not get 2D canvas context.");
}

const atomColorData = [
    [1.0, 0.1, 0.1], // 0: O (Red)
    [0.9, 0.9, 0.9], // 1: H (White)
    [0.2, 0.2, 0.2], // 2: C (Black/Grey)
    [0.1, 0.1, 1.0], // 3: N (Blue)
    [1.0, 1.0, 0.2], // 4: S (Yellow)
    [1.0, 0.5, 0.0], // 5: P (Orange)
    [0.2, 1.0, 0.2], // 6: F (Green)
    [0.5, 0.9, 0.5], // 7: Cl (Light Green)
    [0.6, 0.2, 0.0], // 8: Br (Dark Red)
    [0.4, 0.0, 0.8], // 9: I (Purple)
    [0.2, 1.0, 1.0], // 10: He (Cyan)
];
const atomSymbolData = ["O", "H", "C", "N", "S", "P", "F", "Cl", "Br", "I", "He"];
const atomNames = ["Oxygen", "Hydrogen", "Carbon", "Nitrogen", "Sulfur", "Phosphorus", "Fluorine", "Chlorine", "Bromine", "Iodine", "Helium"];


function getCssColor(type) {
    const color = atomColorData[type];
    if (!color) return 'rgb(255, 0, 255)'; // Magenta for unknown
    const r = Math.floor(color[0] * 255);
    const g = Math.floor(color[1] * 255);
    const b = Math.floor(color[2] * 255);
    return `rgb(${r},${g},${b})`;
}

function drawLegend() {
    ctx.textAlign = 'left';

    const legendX = 15;
    let legendY = 20;
    const lineHeight = 20;
    const swatchRadius = 6;
    const textOffsetX = 25;

    ctx.fillStyle = 'rgba(0, 0, 0, 0.6)';
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.7)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.rect(10, 10, 135, 385);
    ctx.fill();
    ctx.stroke();


    // --- Draw Atom Legend ---
    ctx.fillStyle = 'white';
    ctx.font = 'bold 16px Trebuchet MS';
    ctx.fillText('Legend', legendX + 5, legendY + 10);
    legendY += lineHeight * 1.5;

    ctx.font = '14px Trebuchet MS';
    ctx.fillText('Atoms:', legendX + 5, legendY + 5);
    legendY += lineHeight * 1.5;

    for (let i = 0; i < atomSymbolData.length; i++) {
        // Draw color swatch
        ctx.beginPath();
        ctx.arc(legendX + swatchRadius + 5, legendY - swatchRadius / 2, swatchRadius, 0, 2 * Math.PI);
        ctx.fillStyle = getCssColor(i);
        ctx.fill();

        // Draw text
        ctx.fillStyle = 'white';
        ctx.fillText(`${atomSymbolData[i]} - ${atomNames[i]}`, legendX + textOffsetX, legendY);
        legendY += lineHeight;
    }

    // --- Draw Bond Legend ---
    legendY += lineHeight * 0.5; // Extra space
    ctx.fillText('Bonds:', legendX + 5, legendY);
    legendY += lineHeight;

    // Single Bond
    ctx.beginPath();
    ctx.moveTo(legendX + 5, legendY);
    ctx.lineTo(legendX + 35, legendY);
    ctx.strokeStyle = 'rgba(20, 245, 0, 0.7)';
    ctx.lineWidth = 3;
    ctx.stroke();
    ctx.fillStyle = 'white';
    ctx.fillText('Single', legendX + 45, legendY + 5);
    legendY += lineHeight;

    // Double Bond
    ctx.beginPath();
    ctx.moveTo(legendX + 5, legendY);
    ctx.lineTo(legendX + 35, legendY);
    ctx.strokeStyle = 'rgba(251, 255, 0, 0.85)';
    ctx.lineWidth = 6;
    ctx.stroke();
    ctx.fillText('Double', legendX + 45, legendY + 5);
    legendY += lineHeight;

    // Triple Bond
    ctx.beginPath();
    ctx.moveTo(legendX + 5, legendY);
    ctx.lineTo(legendX + 35, legendY);
    ctx.strokeStyle = 'rgba(255, 0, 0, 0.95)';
    ctx.lineWidth = 9;
    ctx.stroke();
    ctx.fillText('Triple', legendX + 45, legendY + 5);
}


Module.onRuntimeInitialized = () => {
    // --- WebAssembly Function Imports ---
    const initAtoms = Module.cwrap('init_atoms', null, []);
    const simulate = Module.cwrap('simulate', null, ['number']);
    const getAtomsPtr = Module.cwrap('get_atoms', 'number', []);
    const getAtomCount = Module.cwrap('get_atom_count', 'number', []);
    const getAtomSize = Module.cwrap('get_atom_size', 'number', []); // The stride in bytes
    const getAtomTypeOffset = Module.cwrap('get_atom_type_offset', 'number', []);
    const getAtomCovalentRadiusOffset = Module.cwrap('get_atom_covalent_radius_offset', 'number', []);
    const getBondVerticesPtr = Module.cwrap('get_bond_vertices', 'number', []);
    const getVisibleBondCount = Module.cwrap('get_visible_bond_count', 'number', []);
    const spawnAtom = Module.cwrap('spawnAtom', null, ['number', 'number', 'number', 'number', 'number']);

    // --- Simulation Constants from C++ ---
    const SIMULATION_BOX_HALF_WIDTH = 1.0;
    const PICO_TO_NANO = 0.001;
    
    const ATOM_STRUCT_SIZE = getAtomSize();
    const ATOM_TYPE_OFFSET = getAtomTypeOffset();
    const ATOM_RADIUS_OFFSET = getAtomCovalentRadiusOffset();

    // --- Initialize Simulation ---
    initAtoms();

    let atomCount = getAtomCount();
    const atomsPtr = getAtomsPtr();

    // --- Coordinate Transformation Helpers ---
    const simWidth = SIMULATION_BOX_HALF_WIDTH * 2;
    const scale = Math.min(canvas.width, canvas.height) / simWidth;

    const offsetX = (canvas.width - simWidth * scale) / 2;
    const offsetY = (canvas.height - simWidth * scale) / 2;

    function simToCanvasX(simX) {
        return (simX + SIMULATION_BOX_HALF_WIDTH) * scale + offsetX;
    }
    
    function simToCanvasY(simY) {
        return canvas.height - ((simY + SIMULATION_BOX_HALF_WIDTH) * scale + offsetY);
    }
    
    function simToCanvasRadius(simRadius) {
        return simRadius * scale;
    }

    canvas.addEventListener('click', (event) => {
        // Convert canvas coordinates to simulation coordinates
        const rect = canvas.getBoundingClientRect();
        const canvasX = event.clientX - rect.left - offsetX;
        const canvasY = canvas.height - (event.clientY - rect.top - offsetY); // Invert Y for canvas

        // Convert to simulation coordinates
        const simX = (canvasX / scale) - SIMULATION_BOX_HALF_WIDTH;
        const simY = (canvasY / scale) - SIMULATION_BOX_HALF_WIDTH;

        const atomType = Math.floor(Math.random() * atomColorData.length);
        const vx = (Math.random() - 0.5); 
        const vy = (Math.random() - 0.5); 

        spawnAtom(simX, simY, vx, vy, atomType);
    });

    // --- Render Loop ---
    let lastTime = 0;
    function render(time) {
        time *= 0.001; // convert time to seconds
        const deltaTime = time - lastTime;
        lastTime = time;

        simulate(0.016); 

        atomCount = getAtomCount();

        // --- Clear Canvas ---
        ctx.fillStyle = 'rgb(13, 13, 20)'; 
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        const wasmMemory = Module.HEAPU8.buffer;

        // --- RENDER BONDS ---
        const bondCount = getVisibleBondCount();
        if (bondCount > 0) {
            const bondVerticesPtr = getBondVerticesPtr();
            // 5 floats per bond (x1, y1, x2, y2, order)
            const bondData = new Float32Array(wasmMemory, bondVerticesPtr, bondCount * 5);
            
            for (let i = 0; i < bondCount * 5; i += 5) {
                const x1 = simToCanvasX(bondData[i]);
                const y1 = simToCanvasY(bondData[i+1]);
                const x2 = simToCanvasX(bondData[i+2]);
                const y2 = simToCanvasY(bondData[i+3]);
                const order = bondData[i+4];

                switch (order) {
                    case 1:
                        ctx.strokeStyle = 'rgba(20, 245, 0, 0.7)';
                        ctx.lineWidth = 3;
                        break;
                    case 2:
                        ctx.strokeStyle = 'rgba(251, 255, 0, 0.85)';
                        ctx.lineWidth = 6;
                        break;
                    case 3:
                        ctx.strokeStyle = 'rgba(255, 0, 0, 0.95)';
                        ctx.lineWidth = 9;
                        break;
                    default:
                        ctx.strokeStyle = 'rgba(255, 0, 255, 0.7)';
                        ctx.lineWidth = 1;
                        break;
                }
                
                ctx.beginPath();
                ctx.moveTo(x1, y1);
                ctx.lineTo(x2, y2);
                ctx.stroke();
            }
        }

        // --- RENDER ATOMS ---
        for (let i = 0; i < atomCount; i++) {
            const baseOffsetBytes = i * ATOM_STRUCT_SIZE;
            const x = Module.getValue(atomsPtr + baseOffsetBytes + 0, 'float');
            const y = Module.getValue(atomsPtr + baseOffsetBytes + 4, 'float');
            const type = Module.getValue(atomsPtr + baseOffsetBytes + ATOM_TYPE_OFFSET, 'i32');
            const covalentRadiusPM = Module.getValue(atomsPtr + baseOffsetBytes + ATOM_RADIUS_OFFSET, 'i32');
            const simRadius = covalentRadiusPM * PICO_TO_NANO;
            const canvasX = simToCanvasX(x);
            const canvasY = simToCanvasY(y);
            const canvasRadius = simToCanvasRadius(simRadius);

            // Draw the atom circle
            ctx.beginPath();
            ctx.arc(canvasX, canvasY, canvasRadius * 0.93, 0, 2 * Math.PI);
            const atomRgbColor = atomColorData[type];
            ctx.fillStyle = getCssColor(type);
            ctx.fill();

            // Draw Atom Symbol on top
            const symbol = atomSymbolData[type];
            if (symbol) {
                // Determine contrasting text color (black or white) using luminance
                const r = atomRgbColor[0] * 255;
                const g = atomRgbColor[1] * 255;
                const b = atomRgbColor[2] * 255;
                const luminance = (0.299 * r + 0.587 * g + 0.114 * b);
                ctx.fillStyle = luminance > 128 ? 'black' : 'white';

                const fontSize = Math.max(8, canvasRadius); // Minimum size 8px
                ctx.font = `bold ${fontSize}px Arial`;
                ctx.textAlign = 'center';
                ctx.textBaseline = 'middle';
                ctx.fillText(symbol, canvasX, canvasY);
            }
        }
        
        // --- DRAW LEGEND ---
        drawLegend();

        requestAnimationFrame(render);
    }

    requestAnimationFrame(render);
};