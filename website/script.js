/* ===========================
   INTERACTIVE JAVASCRIPT
   =========================== */

// ---- NAV SCROLL EFFECT ----
const nav = document.getElementById('nav');
window.addEventListener('scroll', () => {
    nav.classList.toggle('scrolled', window.scrollY > 60);
});

// ---- ANIMATED COUNTER ----
function animateCounters() {
    document.querySelectorAll('.stat-number').forEach(el => {
        const target = parseInt(el.dataset.target);
        const suffix = el.dataset.suffix || '';
        const duration = 2000;
        const start = performance.now();

        function update(now) {
            const elapsed = now - start;
            const progress = Math.min(elapsed / duration, 1);
            const eased = 1 - Math.pow(1 - progress, 3);
            const current = Math.floor(eased * target);
            el.textContent = current.toLocaleString() + suffix;
            if (progress < 1) requestAnimationFrame(update);
        }
        requestAnimationFrame(update);
    });
}

// ---- INTERSECTION OBSERVER ----
const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.classList.add('visible');

            // Trigger bar animations
            if (entry.target.id === 'importanceChart') {
                entry.target.querySelectorAll('.bar-fill').forEach((bar, i) => {
                    setTimeout(() => bar.classList.add('animated'), i * 100);
                });
            }

            observer.unobserve(entry.target);
        }
    });
}, { threshold: 0.2 });

// Observe elements
document.querySelectorAll('.overview-card, .feature-block, .pipeline-step, .metric-card, .tech-card').forEach(el => {
    el.style.opacity = '0';
    el.style.transform = 'translateY(20px)';
    el.style.transition = 'all 0.6s cubic-bezier(0.25, 0.46, 0.45, 0.94)';
    observer.observe(el);
});

document.querySelectorAll('.visible, .overview-card.visible, .feature-block.visible, .pipeline-step.visible, .metric-card.visible, .tech-card.visible').forEach(el => {
    el.style.opacity = '1';
    el.style.transform = 'translateY(0)';
});

// Re-define visible style via mutation
const visObserver = new MutationObserver((mutations) => {
    mutations.forEach(m => {
        if (m.target.classList.contains('visible')) {
            m.target.style.opacity = '1';
            m.target.style.transform = 'translateY(0)';
        }
    });
});
document.querySelectorAll('.overview-card, .feature-block, .pipeline-step, .metric-card, .tech-card').forEach(el => {
    visObserver.observe(el, { attributes: true, attributeFilter: ['class'] });
});

const chartObserver = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.querySelectorAll('.bar-fill').forEach((bar, i) => {
                setTimeout(() => bar.classList.add('animated'), i * 100);
            });
            chartObserver.unobserve(entry.target);
        }
    });
}, { threshold: 0.3 });

const chart = document.getElementById('importanceChart');
if (chart) chartObserver.observe(chart);

// ---- HERO COUNTER TRIGGER ----
const heroObserver = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            animateCounters();
            heroObserver.unobserve(entry.target);
        }
    });
}, { threshold: 0.5 });

const hero = document.getElementById('hero');
if (hero) heroObserver.observe(hero);

// ---- NETWORK CANVAS ANIMATION ----
const canvas = document.getElementById('networkCanvas');
if (canvas) {
    const ctx = canvas.getContext('2d');
    let width, height, nodes = [], animId;

    function resize() {
        const rect = canvas.parentElement.getBoundingClientRect();
        width = canvas.width = rect.width;
        height = canvas.height = 400;
    }

    function createNodes() {
        nodes = [];
        const count = 40;
        // Some viral nodes (colored differently)
        for (let i = 0; i < count; i++) {
            nodes.push({
                x: Math.random() * width,
                y: Math.random() * height,
                vx: (Math.random() - 0.5) * 0.5,
                vy: (Math.random() - 0.5) * 0.5,
                r: i < 6 ? 6 : 3,
                type: i < 6 ? 'viral' : 'human'
            });
        }
    }

    function draw() {
        ctx.clearRect(0, 0, width, height);

        // Draw edges
        for (let i = 0; i < nodes.length; i++) {
            for (let j = i + 1; j < nodes.length; j++) {
                const dx = nodes[i].x - nodes[j].x;
                const dy = nodes[i].y - nodes[j].y;
                const dist = Math.sqrt(dx * dx + dy * dy);
                if (dist < 100) {
                    const alpha = (1 - dist / 100) * 0.3;
                    // Color based on viral-human interaction
                    if (nodes[i].type !== nodes[j].type) {
                        ctx.strokeStyle = `rgba(0, 212, 255, ${alpha})`;
                    } else {
                        ctx.strokeStyle = `rgba(255, 255, 255, ${alpha * 0.3})`;
                    }
                    ctx.lineWidth = 1;
                    ctx.beginPath();
                    ctx.moveTo(nodes[i].x, nodes[i].y);
                    ctx.lineTo(nodes[j].x, nodes[j].y);
                    ctx.stroke();
                }
            }
        }

        // Draw nodes
        nodes.forEach(node => {
            ctx.beginPath();
            ctx.arc(node.x, node.y, node.r, 0, Math.PI * 2);
            if (node.type === 'viral') {
                ctx.fillStyle = '#00d4ff';
                ctx.shadowColor = '#00d4ff';
                ctx.shadowBlur = 12;
            } else {
                ctx.fillStyle = 'rgba(123, 97, 255, 0.7)';
                ctx.shadowColor = 'transparent';
                ctx.shadowBlur = 0;
            }
            ctx.fill();
            ctx.shadowBlur = 0;

            // Update position
            node.x += node.vx;
            node.y += node.vy;

            // Bounce
            if (node.x < 0 || node.x > width) node.vx *= -1;
            if (node.y < 0 || node.y > height) node.vy *= -1;
        });

        animId = requestAnimationFrame(draw);
    }

    resize();
    createNodes();
    draw();
    window.addEventListener('resize', () => { resize(); createNodes(); });
}

// ---- FLOATING PARTICLES ----
const particleContainer = document.getElementById('particles');
if (particleContainer) {
    const particleCanvas = document.createElement('canvas');
    particleContainer.appendChild(particleCanvas);
    const pCtx = particleCanvas.getContext('2d');

    function resizeParticles() {
        particleCanvas.width = window.innerWidth;
        particleCanvas.height = window.innerHeight;
    }

    resizeParticles();
    window.addEventListener('resize', resizeParticles);

    const particles = [];
    for (let i = 0; i < 50; i++) {
        particles.push({
            x: Math.random() * particleCanvas.width,
            y: Math.random() * particleCanvas.height,
            r: Math.random() * 2 + 0.5,
            vx: (Math.random() - 0.5) * 0.3,
            vy: (Math.random() - 0.5) * 0.3,
            alpha: Math.random() * 0.5 + 0.1
        });
    }

    function drawParticles() {
        pCtx.clearRect(0, 0, particleCanvas.width, particleCanvas.height);
        particles.forEach(p => {
            pCtx.beginPath();
            pCtx.arc(p.x, p.y, p.r, 0, Math.PI * 2);
            pCtx.fillStyle = `rgba(0, 212, 255, ${p.alpha})`;
            pCtx.fill();

            p.x += p.vx;
            p.y += p.vy;

            if (p.x < 0) p.x = particleCanvas.width;
            if (p.x > particleCanvas.width) p.x = 0;
            if (p.y < 0) p.y = particleCanvas.height;
            if (p.y > particleCanvas.height) p.y = 0;
        });
        requestAnimationFrame(drawParticles);
    }
    drawParticles();
}

// ---- SVG GRADIENT FOR RINGS ----
document.querySelectorAll('.metric-ring svg').forEach(svg => {
    const defs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');
    const gradient = document.createElementNS('http://www.w3.org/2000/svg', 'linearGradient');
    gradient.setAttribute('id', 'ringGradient');
    gradient.setAttribute('x1', '0%');
    gradient.setAttribute('y1', '0%');
    gradient.setAttribute('x2', '100%');
    gradient.setAttribute('y2', '0%');

    const stop1 = document.createElementNS('http://www.w3.org/2000/svg', 'stop');
    stop1.setAttribute('offset', '0%');
    stop1.setAttribute('stop-color', '#00d4ff');

    const stop2 = document.createElementNS('http://www.w3.org/2000/svg', 'stop');
    stop2.setAttribute('offset', '100%');
    stop2.setAttribute('stop-color', '#7b61ff');

    gradient.appendChild(stop1);
    gradient.appendChild(stop2);
    defs.appendChild(gradient);
    svg.insertBefore(defs, svg.firstChild);
});
