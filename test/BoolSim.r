library(funcsKishore)
library(BoolNet)

read_boolean_rules <- function(filename) {
  rules_raw <- read_lines(filename) %>%
    .[!str_detect(., "^#|^$")]
  
  rules_df <- tibble(line = rules_raw) %>%
    separate(line, into = c("target", "rule"), sep = " = ", extra = "merge") %>%
    mutate(
      target = str_trim(target),
      rule = str_trim(rule)
    )
  
  return(rules_df)
}

# Convert your rules to BoolNet format
convert_to_boolnet <- function(rules_file) {
  # rules_df should have columns: target, rule
  rules_df <- read_boolean_rules(rules_file)
  boolnet_rules <- rules_df %>%
    mutate(
      # Convert operators to BoolNet syntax
      rule_converted = rule %>%
        str_to_upper() %>%
        str_replace_all("AND", "&") %>%
        str_replace_all("OR", "|") %>%
        str_replace_all("NOT", "!") %>%
        str_replace_all("\\s+", " ")  # Clean up spaces
    ) %>%
    mutate(
      # Format as "target, rule"
      boolnet_line = paste0(target %>% toupper, ", ", rule_converted)
    )
  
  # Write to file
  writeLines(c("targets, factors", boolnet_rules$boolnet_line), str_replace(rules_file, ".txt", ".bnet"))
  
  # cat("BoolNet network written to:", output_file, "\n")
  # return(output_file)
}

plot_boolnet_attractors <- function(attractors, 
                                   node_order = NULL,
                                   colors = c("0" = "blue", "1" = "red"),
                                   title = "Boolean Network Attractors",
                                   attractor_labels = TRUE,
                                   basin_info = T,
                                   text_size = 10,
                                   separator_color = "black",
                                   separator_size = 3) {
  
  # Extract attractor states
  n_attractors <- length(attractors$attractors)
  
  if (n_attractors == 0) {
    stop("No attractors found in the object")
  }
  attractor_data <- plotAttractors(attractors) %>% 
    lapply(function(x) {x %>% t %>% data.frame}) %>% 
    bind_rows %>% 
    mutate(Attractor = rownames(.)) %>% 
    separate(Attractor, into = c("Attractor", "State"), convert = T)
  b <- attractors$stateInfo$attractorAssignment %>% table %>% as.vector
  b <- b/sum(b)
  names(b) <- paste0("Attr", 1:length(b))
  attractor_data$BasinSize <- b[attractor_data$Attractor]
  # Combine all attractors
  combined_data <- attractor_data
  
  # Get node names
  node_names <- colnames(combined_data)
  node_names <- node_names[!node_names %in% c("State", "Attractor", "BasinSize")]
  
  # Apply node ordering if provided
  if (!is.null(node_order)) {
    # Check if all nodes in order exist
    missing_nodes <- setdiff(node_order, node_names)
    if (length(missing_nodes) > 0) {
      warning("Some nodes in node_order not found: ", paste(missing_nodes, collapse = ", "))
    }
    # Use provided order, but include any missing nodes at the end
    extra_nodes <- setdiff(node_names, node_order)
    node_names <- c(intersect(node_order, node_names), extra_nodes)
  }
  
  # Reshape to long format for ggplot
  plot_data <- combined_data %>%
    pivot_longer(
      cols = all_of(node_names),
      names_to = "Node",
      values_to = "Value"
    ) %>%
    mutate(
      Node = factor(Node, levels = rev(node_names)),  # Reverse for top-to-bottom
      Value = factor(Value, levels = c(0, 1)),
      # Create unique state ID across attractors
      GlobalState = paste0("A", Attractor, "_S", State)
    )
  
  # Calculate positions for attractor separators
  state_counts <- combined_data %>%
    group_by(Attractor) %>%
    summarize(n_states = n(), .groups = "drop") %>%
    mutate(
      cumsum_states = cumsum(n_states),
      separator_pos = cumsum_states + 0.5
    )
  
  # Create attractor labels with basin info if requested
  if (attractor_labels) {
    attractor_label_data <- combined_data %>%
      group_by(Attractor) %>%
      summarize(
        n_states = n(),
        mid_state = mean(State),
        BasinSize = first(BasinSize),
        .groups = "drop"
      ) %>%
      mutate(
        cumsum_prev = lag(cumsum(n_states), default = 0),
        label_x = cumsum_prev + (n_states / 2)
      )
    
    if (basin_info && "BasinSize" %in% colnames(attractor_label_data)) {
      attractor_label_data <- attractor_label_data %>%
        mutate(label = paste0("Attractor ", Attractor, 
                            "\n(Basin: ", BasinSize, ")"))
    } else {
      attractor_label_data <- attractor_label_data %>%
        mutate(label = paste0("Attractor ", Attractor))
    }
  }
  
  # Create x-axis positions
  plot_data <- plot_data %>%
    group_by(Attractor) %>%
    mutate(
      states_before = sum(state_counts$n_states[state_counts$Attractor < Attractor[1]]),
      x_position = states_before + State
    ) %>%
    ungroup()
  
  # Create base plot
  p <- ggplot(plot_data %>% 
    mutate(Attractor = paste0("Basin\n", round(BasinSize, 2))), aes(x = x_position, y = Node, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = colors,
      labels = c("0" = "OFF (0)", "1" = "ON (1)"),
      name = "State"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0), add = c(0, 0))) +
    scale_y_discrete(expand = expansion(mult = c(0, 0), add = c(0, 0))) +
    theme_Publication() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = text_size),
      panel.grid = element_blank(),
      legend.position = "top"
    ) +
    facet_grid(
      . ~ Attractor,           # Facet by attractor
      scales = "free_x",       # Allow different x-axis ranges
      space = "free_x"         # KEY: Scale panel width by number of x values
    ) +
    labs(
    #   title = title,
      y = "Nodes"
    )
  
  # Add vertical separators between attractors
#   if (n_attractors > 1) {
#     p <- p +
#       geom_vline(
#         data = state_counts %>% filter(Attractor < n_attractors),
#         aes(xintercept = cumsum_states),
#         color = separator_color,
#         linewidth = separator_size,
#         linetype = "solid"
#       )
#   }
  
  # Add attractor labels
#   if (attractor_labels) {
#     p <- p +
#       geom_text(
#         data = attractor_label_data,
#         aes(x = label_x, y = length(node_names) + 1, label = label),
#         inherit.aes = FALSE,
#         vjust = 0,
#         size = 3.5,
#         fontface = "bold"
#       ) +
#       coord_cartesian(clip = "off") +
#       theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10))
#   }
  
  return(p)
}

# Helper function to get basin information
add_basin_info <- function(attractors) {
  if (is.null(attractors$stateInfo)) {
    message("No basin information available. Run with method='random' or 'exhaustive'")
    return(NULL)
  }
  
  basin_summary <- table(attractors$stateInfo$basinOfAttraction)
  basin_df <- data.frame(
    Attractor = as.integer(names(basin_summary)),
    BasinSize = as.integer(basin_summary)
  )
  
  if (!is.null(attractors$stateInfo$table)) {
    total_states <- nrow(attractors$stateInfo$table)
    basin_df$BasinProportion <- basin_df$BasinSize / total_states
  }
  
  return(basin_df)
}

# Example usage function
boolnet_visualization <- function(rulez, nodeOrder = NULL) {
  # Create example network
  if (!str_detect(rulez, ".bnet$")) {
    convert_to_boolnet(rulez)
    net <- loadNetwork(rulez %>% str_replace(".txt", ".bnet"))
  }
  else {
    net <- loadNetwork(rulez)
  }
  
  nNodes <- length(net$genes)
  # Get attractors with random sampling
  if (nNodes < 14)
    attr <- getAttractors(net, type = "synchronous", method = "exhaustive")
else {
    attr <- getAttractors(net, type = "synchronous", method = "random", 
        startStates = 10000)
}
  
  # Plot with default ordering
  if (is.null(nodeOrder)) {
    p1 <- plot_boolnet_attractors(attr, title = "Default Node Order")
  }
  else {
    p1 <- plot_boolnet_attractors(
        attr, 
        node_order = nodeOrder,
        title = "Custom Node Order",
        attractor_labels = TRUE,
        basin_info = TRUE
    )
  }
  nStates <- attr$attractors %>% 
    sapply(function(x) {
        length(x$involvedStates)
    }) %>% sum
  # Plot with custom ordering
  
  ggsave(str_replace(rulez, ".txt", "_states.png"), width = 3 + 0.3*nStates, height = 2 + 0.2*nNodes, plot = p1)
  ggsave(str_replace(rulez, ".txt", "_states_fixed.png"), width = 10, height = 7, plot = p1)
  # Get basin information
  basin_info <- add_basin_info(attr)
  print("Basin of Attraction Summary:")
  print(basin_info)
  
  return(list(network = net, attractors = attr, basin_info = basin_info))
}

# d <- boolnet_visualization("PSO.txt")
# d <- boolnet_visualization("restriction_switch.txt")
# d <- boolnet_visualization("phase_switch.txt")
# d <- boolnet_visualization("EMT_Switch.txt")
d <- boolnet_visualization("EMT_Switch_Mes.txt")
# d <- boolnet_visualization("Apoptotic_Switch.txt")
# d <- boolnet_visualization("Phase_Switch_NEW.txt")
