/* ----------------------------------------------------------------------------
                | Graph methods
                |---------------------------------------------------------------------------- */

/**
 * Get names of top level keys.
 *
 * @param {object} data
 */
const getSampleNames = (data) =>
    [...Object.keys(data)]

/**
 * Calculates total observation counts per taxa for all samples.
 *
 * @param {object} data
 */
const getSampleCounts = (data) => {
    const counts = {}
    const names = getSampleNames(data)

    const _getSampleCounts = (sample_name, agg, fragment) => {
        Object.entries(fragment).forEach(([key, val]) => {
            agg[key] = {
                [sample_name]: val.count,
                rank: val.rank,
                name: key,
                ...agg[key]
            }
            _getSampleCounts(sample_name, agg, val.children)
        })
    }

    names.map(Sample => _getSampleCounts(Sample, counts, data[Sample]))
    return counts
}

/**
 * Calculates total observation count.
 *
 * @param {object} sample_data
 */
const getTotalCount = (sample_data) =>
    Object.values(sample_data).reduce((sum, I) => I.count + sum, 0)

/**
 * Creates sankey input links list.
 *
 * @param {object} sample_data
 */
const getSankeyEdges = (sample_data) => {
    const _getSankeyEdges = (entries, current) =>
        Object.entries(entries).reduce((arr, [key, val]) => [
            ...arr,
            {
                source: current,
                target: current === key ? `${key}_${val.rank}` : key,
                value: val.count,
                targetRank: val.rank,
                targetRankIdx: ranks.indexOf(val.rank)
            },
            ..._getSankeyEdges(val.children, key)
        ], [])
    return Object.entries(sample_data).map(([k, v]) =>
        _getSankeyEdges(v.children, k)).flat()
}

/**
 * Creates sankey input node list.
 *
 * @param {object} edges.
 */
const getSankeyNodes = (edges) => {
    const unique = new Set(
        edges.reduce((arr, Link) =>
            [Link.source, Link.target, ...arr], [])
    )
    return [...unique].map(Item => ({ name: Item }))
}

/**
 * Creates a sankey graph generator.
 *
 * @param {number} width
 * @param {number} height
 */
const getSankeyGenerator = (width, height) =>
    d3.sankey()
        .nodeId(d => d.name)
        // .nodeSort(() => true) # Disabled for now
        .nodeAlign(d3['sankeyCenter'])
        .nodeWidth(30)
        .nodePadding(30)
        .extent([[0, 5], [width, height - 5]])

/**
 * Creates a sankey graph.
 *
 * @param {string} sample_name
 * @param {object} sample_data
 * @param {object} generator
 * @param {number} cutoff
 * @param {number} total
 */
const getSankeyGraph = (sample_name, sample_data, generator, cutoff, total) => {
    const _cutoffVal = total / 100 * cutoff

    const _edges = getSankeyEdges(sample_data)
    const _filteredEdges = _edges.filter(Link => Link.value >= _cutoffVal)
    const _nodes = getSankeyNodes(_filteredEdges)

    return generator({
        nodes: _nodes.map(d => Object.assign({}, d)),
        links: _filteredEdges.map(d => Object.assign({}, d))
    });
}

/* ----------------------------------------------------------------------------
| Render chart
|---------------------------------------------------------------------------- */

/**
 * Renders a selection dropdown.
 *
 * @param {string} select_id
 * @param {array} sample_names
 * @param {string} symbol
 */
const renderSelect = (select_id, sample_names, symbol = null) => {
    const dropdown = d3.select(select_id)
    dropdown
        .selectAll('option')
        .data(sample_names)
        .enter()
        .append("option")
        .attr("value", (d) => d)
        .text((d) => {
            return symbol ? d + symbol : d;
        })
    return dropdown
}

/**
 * Renders a sankey node tooltip
 */
const renderToolTop = () =>
    d3.select("#tooltip")
        .style("position", "absolute")
        .style("z-index", "10")
        .style("visibility", "hidden")
        .style("background", "white");

/**
 * Sets the text style for a given selection.
 *
 * @param {object} enter
 * @param {string} colour
 */
const setTextStyles = (enter, colour = "#555") => {
    enter.style("font-family", "monospace")
        .style("text-transform", "lowercase")
        .style("letter-spacing", "0.06em")
        .attr("fill", "#555")
}

/**
 * Updates a node tooltip on mouse events
 *
 * @param {object} enter
 * @param {object} colourScale
 * @param {number} total
 * @param {object} counts
 */
const updateToolTip = (enter, colourScale, total, counts) => {
    const toolTip = d3.select("#tooltip")
    enter.on("mouseover", (event, d) => {
        const toolTipData = [
            `Name: ${d.name}`,
            `Rank: ${counts[d.name].rank}`,
            `Count: ${d.value}`,
            `Percentage: ${(100 / total * d.value).toFixed(2)}%`
        ]
        toolTip
            .select("text")
            .remove()
        const toolTipText = toolTip
            .style("top", d.y + "px")
            .style("left", d.x + "px")
            .style("padding", "5px")
            .style("border-radius", "4px")
            .style("visibility", "visible")
            .style("background-color", "#fafafa")
            .style("border", `2px solid ${colourScale(d.value)}`)
            .style("box-shadow", "5px 12px 20px rgb(36 37 38 / 13%)")
            .append('text')
            .attr("x", 0)
            .attr("y", 0)
            .call(enter => setTextStyles(enter))
        toolTipText
            .selectAll("tspan")
            .data(toolTipData)
            .join("tspan")
            .attr("x", 0)
            .attr("y", 0)
            .attr("dy", (d, i) => `${1.2 * i}em`)
            .style("display", "block")
            .attr("fill-opacity", 0.7)
            .text(d => d)
    })
        .on("mousemove", (event, d) => {
            toolTip
                .style("top", event.layerY + 10 + "px")
                .style("left", event.layerX + 10 + "px")
        })
        .on("mouseout", (d) => {
            return toolTip.style("visibility", "hidden")
        });
}

/**
 * Updates the sankey plot.
 *
 * @param {object} svg
 * @param {object} graph
 * @param {number} depth
 * @param {object} colourScale
 * @param {number} total
 * @param {object} counts
 */
const updateVisualisation = (svg, graph, depth, colourScale, total, counts) => {
    const { nodes, links } = graph

    const _filteredNodes = nodes.filter(Node =>
        ranks.indexOf(counts[Node.name]?.rank || 0) <= depth
        && Node.name !== 'Unknown')

    // Temporary extra jank: Remove the node and link for unknown species
    const _filteredLinks = links.filter(Link =>
        Link.targetRankIdx <= depth && Link.target.name !== 'Unknown')

    const t = svg.transition().duration(750);

    const nodelist = svg.select("#nodes")
        .selectChildren("g")
        .data(_filteredNodes, (d) => d.name)
        .join(
            enter => {
                const container = enter.append("g")

                // Add node rect
                container
                    .append("rect")
                    .attr("x", d => d.x0 + 1)
                    .attr("y", d => d.y0)
                    .attr("height", 0)
                    .attr("width", d => d.x1 - d.x0 - 2)
                    .attr("fill", d => colourScale(d.value))
                    .style("cursor", "pointer")
                    .attr("pointer-events", "all")
                    .call(enter => enter.transition(t)
                        .attr("height", d => d.y1 - d.y0))
                    .call(enter => updateToolTip(enter, colourScale, total, counts))

                // Add node accessibility title
                container
                    .append("title")
                    .text(d => `${d.name}\n${d.value.toLocaleString()}`)

                // Add node label
                container
                    .append("text")
                    .attr("x", d => d.x1 + 6)
                    .attr("y", d => (d.y1 + d.y0) / 2)
                    .attr("height", d => d.y1 - d.y0)
                    .attr("width", d => d.x1 - d.x0 - 2)
                    .attr("dy", "0.35em")
                    .attr("text-anchor", "start")
                    .text(d => d.name)
                    .call(enter => setTextStyles(enter))
                    .append("tspan")
                    .attr("dy", "1.2em")
                    .attr("x", d => d.x1 + 6)
                    .attr("fill-opacity", 0.7)
                    .text(d => `${d.value.toLocaleString()}`);
            },
            update => {
                update.select("rect")
                    .attr("fill", d => colourScale(d.value))
                    .call(enter => enter.transition(t)
                        .attr("x", d => d.x0 + 1)
                        .attr("y", d => d.y0)
                        .attr("height", d => d.y1 - d.y0)
                        .attr("width", d => d.x1 - d.x0 - 2))

                // Add node accessibility title
                update.select("title")
                    .text(d => `${d.name}\n${d.value.toLocaleString()}`)

                update.select("text")
                    .attr("x", d => d.x1 + 6)
                    .attr("y", d => (d.y1 + d.y0) / 2)
                    .attr("height", d => d.y1 - d.y0)
                    .attr("width", d => d.x1 - d.x0 - 2)
                    .attr("dy", "0.35em")
                    .attr("text-anchor", "start")
                    .text(d => d.name)
                    .append("tspan")
                    .attr("dy", "1.2em")
                    .attr("x", d => d.x1 + 6)
                    .attr("fill-opacity", 0.7)
                    .text(d => `${d.value.toLocaleString()}`);
            },
            exit => {
                exit.remove()
            }
        )

    const linkList = svg.select("#links")
        .selectChildren("g")
        .data(_filteredLinks, (d) => `${d.source.name}-${d.target.name}`)
        .join(
            enter => {
                const container = enter.append("g")
                    .style("mix-blend-mode", "multiply")

                container
                    .append("linearGradient")
                    .attr("id", d => `${d.source.name}-${normaliseName(d.target.name)}-grad`)
                    .attr("gradientUnits", "userSpaceOnUse")
                    .attr("x1", d => d.source.x1)
                    .attr("x2", d => d.target.x0)
                    .call(gradient => gradient.append("stop")
                        .attr("offset", "0%")
                        .attr("stop-color", (d) => colourScale(d.source.value)))
                    .call(gradient => gradient.append("stop")
                        .attr("offset", "100%")
                        .attr("stop-color", (d) => colourScale(d.target.value)));

                container
                    .append("path")
                    .attr("d", d3.sankeyLinkHorizontal())
                    .attr("stroke", (d) => `url(#${d.source.name}-${normaliseName(d.target.name)}-grad`)
                    .attr("stroke-opacity", 0.1)
                    .call(enter => enter.transition(t)
                        .attr("stroke-width", d => Math.max(1, d.width)))
            },
            update => {
                update
                    .attr("stroke", d => colourScale(d.target.value))
                    .style("mix-blend-mode", "multiply")

                update
                    .select("path")
                    .attr("d", d3.sankeyLinkHorizontal())
                    .call(enter => enter.transition(t)
                        .attr("stroke-width", d => Math.max(1, d.width)))
            },
            exit => {
                exit.remove()
            }
        )

    return svg
}

/**
 * Renders the initial plot.
 *
 * @param {object} svg
 * @param {object} graph
 * @param {number} depth
 * @param {object} colourScale
 * @param {number} total
 * @param {object} counts
 */
const renderVisualisation = (svg, graph, depth, colourScale, total, counts) => {
    const nodelist = svg.append("g")
        .attr("id", "nodes")
    const linklist = svg.append("g")
        .attr("id", "links")
        .attr("fill", "none")

    updateVisualisation(svg, graph, depth, colourScale, total, counts)
    return svg
}

/* ----------------------------------------------------------------------------
| Chart reactivity
|---------------------------------------------------------------------------- */

/**
 * Handles updating the plot when a setting is changed
 *
 * @param {object} counts
 */
const handlePlotSelectChange = (counts) => {
    const rank = d3.select("#rank-select").property('value');
    const cutoff = d3.select("#cutoff-select").property('value');
    const sample = d3.select("#sample-select").property('value');
    const sample_data = parsed[sample]
    const _total = getTotalCount(sample_data)
    const _graph = getSankeyGraph(sample, parsed[sample], generator, cutoff, _total)
    const colourScale = d3.scaleQuantize().domain([0, _total]).range(colours);
    setStateGraph(_graph)
    updateVisualisation(svg, _graph, ranks.indexOf(rank), colourScale, _total, counts);
}

/**
 * Handles zooming on the sankey plot
 *
 * @param {object} e
 */
const handleZoom = (e) => {
    d3.select('svg')
        .selectChildren('g')
        .attr('transform', e.transform)
    const fontSize = Math.min(16 / e.transform.k, 12)
    d3.selectAll('#nodes text')
        .style("font-size", `${fontSize}px`)
}

/* ----------------------------------------------------------------------------
| Render table
|---------------------------------------------------------------------------- */

/**
 * Renders the table rank selector
 *
 * @param {string} _id
 * @param {array} ranks
 */
const renderTableRankSelect = (_id, ranks) => {
    const rankSelect = d3.select('.dataTable-top')
        .insert('div', ":first-child")
        .classed("dataTable-rank", true)
        .append('label')
        .text('Select rank')

    rankSelect
        .insert('select')
        .attr('id', _id)

    return renderSelect(`#${_id}`, ranks)
}

/**
 * Renders the initial table
 *
 * @param {object} counts
 * @param {array} samples
 * @param {string} rank
 */
const renderTable = (counts, samples, rank) => {
    const table = d3.select('#table')

    const headerRows = ['Taxon', 'Rank', 'Total', ...samples]
    const thead = table
        .select("thead tr")
        .selectAll("th")
        .data(headerRows)
        .join('th')
        .text(d => d)

    const bodyRows = Object.values(counts).filter(Count => Count.rank === rank)

    const trows = table
        .select("tbody")
        .selectAll("tr")
        .data(bodyRows)
        .join('tr')
        .selectAll('td')
        .data((d) => {
            const sampleCounts = samples.map(Sample => d[Sample] || 0)
            const total = sampleCounts.reduce((sum, I) => I + sum, 0)
            return [d.name, d.rank, total, ...sampleCounts]
        })
        .join('td')
        .text((d) => d)
}

/**
 * Updates the counts table on rank select change
 *
 * @param {object} datatable
 * @param {object} counts
 * @param {array} samples
 */
const handleTableSelectChange = (datatable, counts, samples) => {
    const table = d3.select('#table')
    const rank = d3.select("#table-rank-select").property('value');

    datatable.rows().remove(
        [...datatable.data.map(R => R.dataIndex)]
    )

    const bodyRows = Object.values(counts)
        .filter(Count => Count.rank === rank)
        .map(Row => {
            const sampleCounts = samples.map(Sample => (Row[Sample] || 0))
            const total = sampleCounts.reduce((sum, I) => I + sum, 0).toString()
            return [
                Row.name,
                Row.rank,
                total,
                ...sampleCounts.map(Count => Count.toString())
            ]
        })

    datatable.rows().add(bodyRows)
    datatable.refresh();
}

/* ----------------------------------------------------------------------------
                | Initialisation
                |---------------------------------------------------------------------------- */

const svg = d3.select('#sankey-plot')
    .insert('svg')
    .attr("width", width)
    .attr("height", height)
    .style("background", "#fff")

// FYI CONSIDER THIS GLOBAL
const state = {}
const setStateGraph = (graph) => {
    state.graph = graph
}

// Prepare data
const parsed = parseData(getData())
const generator = getSankeyGenerator(width, height)

// Initialise samples
const names = getSampleNames(parsed)
const default_sample = names[0]
const sample_data = parsed[default_sample]
const sample_counts = getSampleCounts(parsed)
const sampleSelect = renderSelect('#sample-select', names, 0)
sampleSelect.property('value', default_sample);

// Initialise ranks
const default_rank = "order"
const rankSelect = renderSelect('#rank-select', ranks)
rankSelect.property('value', default_rank);

// Initialise cutoff
const default_cutoff = 1
const cutoffSelect = renderSelect('#cutoff-select', cutoffs, '%')
cutoffSelect.property('value', `${default_cutoff}`);

// Render default sample
const total = getTotalCount(sample_data)
const graph = getSankeyGraph(
    default_sample, sample_data, generator, default_cutoff, total)
const colourScale = d3.scaleQuantize().domain([0, total]).range(colours);
renderVisualisation(svg, graph, ranks.indexOf(default_rank), colourScale, total, sample_counts)
setStateGraph(graph)

// Initialise select reactivity
d3.select("#sample-select").on("change", () => handlePlotSelectChange(sample_counts));
d3.select("#rank-select").on("change", () => handlePlotSelectChange(sample_counts));
d3.select("#cutoff-select").on("change", () => handlePlotSelectChange(sample_counts));

// Initialise zoom reactivity
const zoom = d3.zoom().on('zoom', handleZoom);
svg.call(zoom);

// Initialise tooltip interactivity
const toolTip = renderToolTop()

// Initialise manual zoom interactivity
function zoomIn() {
    d3.select('svg')
        .transition()
        .call(zoom.scaleBy, 2);
}

function zoomOut() {
    d3.select('svg')
        .transition()
        .call(zoom.scaleBy, 0.5);
}

function resetZoom() {
    d3.select('svg')
        .transition()
        .call(zoom.scaleTo, 1);
}

// Render table
renderTable(sample_counts, names, default_rank)
const dataTable = new simpleDatatables.DataTable("#table", {
    searchable: true,
    fixedHeight: true,
    columns: [
        { select: 2, sort: "desc" },
    ]
})

// Initialise chart interactivity
const tableRankSelect = renderTableRankSelect('table-rank-select', ranks)
tableRankSelect.property('value', default_rank);
d3.select("#table-rank-select").on("change", () => handleTableSelectChange(
    dataTable, sample_counts, names));