// ---- PREDICTION ----
async function runPrediction() {
    const viral = document.getElementById('viralInput').value.trim();
    const human = document.getElementById('humanInput').value.trim();
    const btn = document.getElementById('predictBtn');
    const btnText = document.getElementById('btnText');
    const spinner = document.getElementById('btnSpinner');
    const errorEl = document.getElementById('errorMsg');
    const resultEl = document.getElementById('resultBox');

    if (!viral || !human) {
        errorEl.textContent = 'Please enter both protein IDs.';
        errorEl.classList.remove('hidden');
        resultEl.classList.add('hidden');
        return;
    }

    // Loading state
    btn.disabled = true;
    btnText.textContent = 'Predicting...';
    spinner.classList.remove('hidden');
    errorEl.classList.add('hidden');
    resultEl.classList.add('hidden');

    try {
        const resp = await fetch('/predict', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ viral_id: viral, human_id: human })
        });

        const data = await resp.json();

        if (!resp.ok) {
            throw new Error(data.error || 'Prediction failed');
        }

        // Show result
        const badge = document.getElementById('resultBadge');
        badge.textContent = data.prediction;
        badge.className = 'badge ' + (data.prediction === 'INTERACTING' ? 'interacting' : 'non-interacting');

        document.getElementById('probValue').textContent = data.probability.toFixed(4);
        document.getElementById('confValue').textContent = data.confidence.toFixed(1);

        const bar = document.getElementById('resultBar');
        bar.style.width = (data.probability * 100) + '%';

        // Viral details
        document.getElementById('vId').textContent = data.viral.id;
        document.getElementById('vName').textContent = data.viral.name;
        document.getElementById('vLen').textContent = data.viral.seq_length + ' aa';
        document.getElementById('vMw').textContent = data.viral.mw.toLocaleString() + ' Da';
        document.getElementById('vGravy').textContent = data.viral.gravy;

        // Human details
        document.getElementById('hId').textContent = data.human.id;
        document.getElementById('hName').textContent = data.human.name;
        document.getElementById('hLen').textContent = data.human.seq_length + ' aa';
        document.getElementById('hMw').textContent = data.human.mw.toLocaleString() + ' Da';
        document.getElementById('hGravy').textContent = data.human.gravy;

        resultEl.classList.remove('hidden');

    } catch (err) {
        errorEl.textContent = err.message;
        errorEl.classList.remove('hidden');
    } finally {
        btn.disabled = false;
        btnText.textContent = 'Predict';
        spinner.classList.add('hidden');
    }
}

// Enter key to submit
document.getElementById('humanInput').addEventListener('keypress', e => {
    if (e.key === 'Enter') runPrediction();
});
document.getElementById('viralInput').addEventListener('keypress', e => {
    if (e.key === 'Enter') document.getElementById('humanInput').focus();
});

// ---- NETWORK GRAPH ----
const canvas = document.getElementById('networkCanvas');
const ctx = canvas.getContext('2d');
let nodes = [], edges = [];
let width, height;
let dragging = null, offsetX = 0, offsetY = 0;
let hovered = null;

function initNetwork() {
    const rect = canvas.parentElement.getBoundingClientRect();
    width = canvas.width = rect.width;
    height = canvas.height = 800;

    const data = NETWORK_DATA;
    if (!data.nodes.length) return;

    // Position nodes - viral in center ring, human around
    const viralNodes = data.nodes.filter(n => n.type === 'viral');
    const humanNodes = data.nodes.filter(n => n.type === 'human');

    const cx = width / 2, cy = height / 2;

    viralNodes.forEach((n, i) => {
        const angle = (i / viralNodes.length) * Math.PI * 2;
        const r = 180;
        nodes.push({
            ...n, x: cx + Math.cos(angle) * r, y: cy + Math.sin(angle) * r,
            vx: 0, vy: 0, radius: 18
        });
    });

    humanNodes.forEach((n, i) => {
        const angle = (i / humanNodes.length) * Math.PI * 2 + 0.5;
        const r = 250 + Math.random() * 150;
        nodes.push({
            ...n, x: cx + Math.cos(angle) * r, y: cy + Math.sin(angle) * r,
            vx: 0, vy: 0, radius: 8
        });
    });

    edges = data.edges.map(e => ({
        source: nodes.find(n => n.id === e.source),
        target: nodes.find(n => n.id === e.target)
    })).filter(e => e.source && e.target);

    // Simple force simulation
    for (let iter = 0; iter < 120; iter++) {
        // Repulsion
        for (let i = 0; i < nodes.length; i++) {
            for (let j = i + 1; j < nodes.length; j++) {
                const dx = nodes[j].x - nodes[i].x;
                const dy = nodes[j].y - nodes[i].y;
                const dist = Math.max(1, Math.sqrt(dx * dx + dy * dy));
                const force = 1200 / (dist * dist);
                nodes[i].vx -= (dx / dist) * force;
                nodes[i].vy -= (dy / dist) * force;
                nodes[j].vx += (dx / dist) * force;
                nodes[j].vy += (dy / dist) * force;
            }
        }
        // Attraction (edges)
        edges.forEach(e => {
            const dx = e.target.x - e.source.x;
            const dy = e.target.y - e.source.y;
            const dist = Math.sqrt(dx * dx + dy * dy);
            const force = (dist - 120) * 0.01;
            e.source.vx += (dx / dist) * force;
            e.source.vy += (dy / dist) * force;
            e.target.vx -= (dx / dist) * force;
            e.target.vy -= (dy / dist) * force;
        });
        // Center gravity
        nodes.forEach(n => {
            n.vx += (cx - n.x) * 0.005;
            n.vy += (cy - n.y) * 0.005;
            n.x += n.vx * 0.3;
            n.y += n.vy * 0.3;
            n.vx *= 0.8;
            n.vy *= 0.8;
            // Keep in bounds
            n.x = Math.max(20, Math.min(width - 20, n.x));
            n.y = Math.max(20, Math.min(height - 20, n.y));
        });
    }

    drawNetwork();
}

function drawNetwork() {
    ctx.clearRect(0, 0, width, height);

    // Draw edges
    edges.forEach(e => {
        const isHighlighted = hovered && (hovered.id === e.source.id || hovered.id === e.target.id);
        ctx.beginPath();
        ctx.moveTo(e.source.x, e.source.y);
        ctx.lineTo(e.target.x, e.target.y);
        ctx.strokeStyle = isHighlighted ? 'rgba(6,182,212,0.5)' : 'rgba(255,255,255,0.06)';
        ctx.lineWidth = isHighlighted ? 1.5 : 0.5;
        ctx.stroke();
    });

    // Draw nodes
    nodes.forEach(n => {
        const isHov = hovered && hovered.id === n.id;
        const isConnected = hovered && edges.some(e =>
            (e.source.id === hovered.id && e.target.id === n.id) ||
            (e.target.id === hovered.id && e.source.id === n.id)
        );

        ctx.beginPath();
        ctx.arc(n.x, n.y, isHov ? n.radius + 3 : n.radius, 0, Math.PI * 2);

        if (n.type === 'viral') {
            ctx.fillStyle = isHov ? '#22d3ee' : '#06b6d4';
            ctx.shadowColor = '#06b6d4';
            ctx.shadowBlur = isHov ? 20 : 12;
        } else {
            ctx.fillStyle = isConnected ? '#a78bfa' : (isHov ? '#a78bfa' : 'rgba(139,92,246,0.6)');
            ctx.shadowColor = '#a78bfa';
            ctx.shadowBlur = isConnected || isHov ? 10 : 0;
        }
        ctx.fill();
        ctx.shadowBlur = 0;

        // Labels for viral proteins
        if (n.type === 'viral') {
            ctx.font = '700 12px Inter, sans-serif';
            ctx.textAlign = 'center';
            ctx.fillStyle = '#ffffff';
            ctx.fillText(n.label, n.x, n.y - n.radius - 8);
        }

        // Labels for human proteins (show always, highlight on hover)
        if (n.type === 'human') {
            ctx.font = isHov ? '600 11px Inter, sans-serif' : '500 10px Inter, sans-serif';
            ctx.textAlign = 'center';
            ctx.fillStyle = isHov || isConnected ? '#e5e7eb' : '#9ca3af';
            ctx.fillText(n.id, n.x, n.y - n.radius - 6);
        }
    });

    requestAnimationFrame(drawNetwork);
}

// Mouse interactions for network
canvas.addEventListener('mousemove', e => {
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;

    if (dragging) {
        dragging.x = mx + offsetX;
        dragging.y = my + offsetY;
        return;
    }

    let found = null;
    for (const n of nodes) {
        const dx = mx - n.x, dy = my - n.y;
        if (Math.sqrt(dx * dx + dy * dy) < n.radius + 5) {
            found = n;
            break;
        }
    }
    hovered = found;
    canvas.style.cursor = found ? 'pointer' : 'grab';

    // Update info panel
    const info = document.getElementById('nodeInfo');
    if (found) {
        document.getElementById('nodeInfoTitle').textContent = found.type === 'viral'
            ? `🦠 ${found.label}` : `🧬 ${found.id}`;
        const connCount = edges.filter(e => e.source.id === found.id || e.target.id === found.id).length;
        document.getElementById('nodeInfoDetail').textContent =
            `${found.type === 'viral' ? 'Viral' : 'Human'} protein · ${connCount} interactions`;
        info.classList.remove('hidden');
    } else {
        info.classList.add('hidden');
    }
});

canvas.addEventListener('mousedown', e => {
    if (hovered) {
        dragging = hovered;
        const rect = canvas.getBoundingClientRect();
        offsetX = hovered.x - (e.clientX - rect.left);
        offsetY = hovered.y - (e.clientY - rect.top);
    }
});

canvas.addEventListener('mouseup', () => { dragging = null; });
canvas.addEventListener('mouseleave', () => { dragging = null; hovered = null; });

// Init
window.addEventListener('load', () => {
    initNetwork();
});
window.addEventListener('resize', () => {
    nodes = []; edges = [];
    initNetwork();
});
