(function () {
  "use strict";

  function currentPath() {
    try {
      return new URL(window.location.href).pathname;
    } catch (e) {
      return window.location.pathname || "/";
    }
  }

  function markActivePage() {
    var path = currentPath().replace(/\/+$/, "") || "/";
    document.querySelectorAll(".topnav a").forEach(function (link) {
      var href = "";
      try {
        href = new URL(link.getAttribute("href"), window.location.href).pathname;
      } catch (e) {
        return;
      }
      href = href.replace(/\/+$/, "") || "/";
      if (href === path) {
        link.classList.add("active");
      }
    });
  }

  function updateTaskBadge() {
    var link = document.querySelector('.topnav a[href="/tasks"]');
    if (!link) return;
    var badge = link.querySelector(".nav-count");
    if (!badge) {
      badge = document.createElement("span");
      badge.className = "nav-count";
      link.appendChild(badge);
    }
    fetch("/tasks/count", {headers: {"Accept": "application/json"}})
      .then(function (response) {
        if (!response.ok) throw new Error("bad status");
        return response.json();
      })
      .then(function (data) {
        var count = Number(data && data.count);
        if (!count || count < 0) {
          badge.hidden = true;
          badge.textContent = "";
          return;
        }
        badge.textContent = String(count);
        badge.hidden = false;
      })
      .catch(function () {
        badge.hidden = true;
      });
  }

  markActivePage();
  updateTaskBadge();
  window.setInterval(updateTaskBadge, 4000);
})();
