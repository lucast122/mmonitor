function updateActiveNavLink() {
     var currentPath = window.location.pathname;
     document.querySelectorAll('.nav-link').forEach(function (link) {
          if (currentPath.includes(link.getAttribute('href'))) {
               link.classList.add('active');
          } else {
               link.classList.remove('active');
          }
     });
}

document.addEventListener('DOMContentLoaded', updateActiveNavLink);
window.addEventListener('popstate', updateActiveNavLink); 