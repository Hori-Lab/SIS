// Custom JavaScript for interactive toctree

document.addEventListener('DOMContentLoaded', function() {
    // Add click handlers to toctree items
    const toctreeItems = document.querySelectorAll('.toctree-wrapper li');
    
    toctreeItems.forEach(function(item) {
        const link = item.querySelector('a');
        const subList = item.querySelector('ul');
        
        if (subList) {
            // Add click handler for expand/collapse
            link.addEventListener('click', function(e) {
                e.preventDefault();
                item.classList.toggle('expanded');
                
                // Store state in localStorage
                const isExpanded = item.classList.contains('expanded');
                localStorage.setItem('toctree-' + link.href, isExpanded);
            });
            
            // Restore state from localStorage
            const savedState = localStorage.getItem('toctree-' + link.href);
            if (savedState === 'true') {
                item.classList.add('expanded');
            }
        }
    });
    
    // Add keyboard navigation
    document.addEventListener('keydown', function(e) {
        if (e.key === 'Enter' || e.key === ' ') {
            const focusedItem = document.querySelector('.toctree-wrapper li:focus');
            if (focusedItem) {
                const link = focusedItem.querySelector('a');
                if (link) {
                    link.click();
                }
            }
        }
    });
}); 